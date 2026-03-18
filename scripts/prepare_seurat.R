#!/usr/bin/env Rscript
# scripts/prepare_seurat.R

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
  library(SeuratDisk)
  library(jsonlite)
  library(png)
  library(arrow)
})

option_list <- list(
  make_option(c("--data_dir"), type="character", help="Directory containing input files"),
  make_option(c("--filename"), type="character", default="filtered_feature_bc_matrix.h5", help="HDF5 filename"),
  make_option(c("--out"), type="character", help="Output .rds path"),
  make_option(c("--platform"), type="character", default="Visium", help="Platform: Visium, VisiumHD, SlideSeq, StereoSeq"),
  make_option(c("--bin.size"), type="numeric", default=8, help="Bin size for VisiumHD: 8 or 16")
)

opt <- parse_args(OptionParser(option_list=option_list))

data_dir <- opt$data_dir
filename <- opt$filename
outfile  <- opt$out
platform <- opt$platform
bin_size <- opt$bin.size

if (!dir.exists(data_dir)) stop("data_dir not found: ", data_dir)

# ==============================================================================
# 1. SMART PATH DETECTION
# ==============================================================================
full_h5_path <- file.path(data_dir, filename)

if (platform == "VisiumHD" && !file.exists(full_h5_path)) {
  message("VisiumHD: Checking specified bin size (", bin_size, "um) for H5...")
  bin_string <- sprintf("%03dum", bin_size)
  other_bin <- if(bin_size == 8) "016um" else "008um"
  
  candidates <- c(
    file.path("binned_outputs", paste0("square_", bin_string), "filtered_feature_bc_matrix.h5"),
    file.path("binned_outputs", paste0("square_", other_bin), "filtered_feature_bc_matrix.h5")
  )
  
  found <- FALSE
  for (cand in candidates) {
    if (file.exists(file.path(data_dir, cand))) {
      filename <- cand
      found <- TRUE
      break
    }
  }
  if (!found) {
    all_h5 <- list.files(data_dir, pattern = "filtered_feature_bc_matrix\\.h5$", recursive = TRUE)
    if (length(all_h5) > 0) filename <- all_h5[1]
  }
}
full_h5_path <- file.path(data_dir, filename)
#if (!file.exists(full_h5_path)) stop("Matrix H5 file not found: ", full_h5_path)

message(Sys.time(), " - Detected Matrix: ", full_h5_path)

# ==============================================================================
# 2. ROBUST OBJECT CREATION
# ==============================================================================

if (platform %in% c("Visium", "VisiumHD")) {
  
  message(" - Loading Counts...")
  counts <- Read10X_h5(full_h5_path)
  # Ensure the assay is named "Spatial" from the start
  seu <- CreateSeuratObject(counts = counts, assay = "Spatial")
  
  message(" - Loading Spatial Image/Coords...")
  img_search_paths <- unique(c(
    file.path(data_dir, "spatial"),
    file.path(dirname(full_h5_path), "spatial"),
    dirname(full_h5_path)
  ))
  
  spatial_dir <- NULL
  for (p in img_search_paths) {
    if (file.exists(file.path(p, "scalefactors_json.json")) || 
        length(list.files(p, pattern = "tissue_positions")) > 0) {
      spatial_dir <- p
      break
    }
  }
  
  if (!is.null(spatial_dir)) {
    message("    Found spatial info at: ", spatial_dir)
    
    coord_file <- list.files(spatial_dir, pattern = "tissue_positions", full.names = TRUE)[1]
    if (is.na(coord_file)) {
      coord_file <- list.files(dirname(spatial_dir), pattern = "tissue_positions", full.names = TRUE)[1]
    }
    
    if (!is.na(coord_file)) {
      if (grepl("\\.parquet$", coord_file)) {
        positions <- arrow::read_parquet(coord_file)
        positions <- as.data.frame(positions)
        rownames(positions) <- as.character(positions$barcode)
      } else {
        positions <- read.csv(coord_file, header = FALSE, row.names = 1)
        colnames(positions) <- c("tissue", "row", "col", "imagerow", "imagecol")
      }
      
      # Standardize VisiumHD column names to match legacy VisiumV1 format
      if ("pxl_row_in_fullres" %in% colnames(positions)) {
        message("    Reformatting VisiumHD coordinates to standard format...")
        colnames(positions)[colnames(positions) == "in_tissue"] <- "tissue"
        colnames(positions)[colnames(positions) == "array_row"] <- "row"
        colnames(positions)[colnames(positions) == "array_col"] <- "col"
        colnames(positions)[colnames(positions) == "pxl_row_in_fullres"] <- "imagerow"
        colnames(positions)[colnames(positions) == "pxl_col_in_fullres"] <- "imagecol"
        positions <- positions[, c("tissue", "row", "col", "imagerow", "imagecol")]
      }
      
      common_barcodes <- intersect(colnames(seu), rownames(positions))
      if (length(common_barcodes) > 0) {
        positions <- positions[common_barcodes, ]
        
        sf_path <- file.path(spatial_dir, "scalefactors_json.json")
        if (file.exists(sf_path)) {
          sf <- fromJSON(sf_path)
          
          # Safe retrieval of scale factors
          lowres_sf <- if (!is.null(sf$tissue_lowres_scalef)) sf$tissue_lowres_scalef else 1
          
          lowres_path <- file.path(spatial_dir, "tissue_lowres_image.png")
          img_data <- if (file.exists(lowres_path)) png::readPNG(lowres_path) else matrix(1, 1, 1)
          
          img <- new(Class = "VisiumV1", 
                     image = img_data, 
                     scale.factors = scalefactors(
                       spot = sf$spot_diameter_fullres,
                       fiducial = sf$fiducial_diameter_fullres,
                       hires = sf$tissue_hires_scalef,
                       lowres = lowres_sf
                     ),
                     coordinates = positions,
                     spot.radius = (sf$spot_diameter_fullres * lowres_sf))
          
          # --- THE FIX ---
          # Explicitly link the image to the "Spatial" assay
          img@assay <- "Spatial"
          img@key <- "slice1_"
          
          # Now this assignment will work without the "length zero" error
          seu[["slice1"]] <- img
        }
      }
    }
  }
} else if (platform == "Slideseq") {
  mtx_file <- list.files(data_dir, pattern = "\\.mtx$", full.names = TRUE)
  csv_file <- list.files(data_dir, pattern = "BeadLocations.*\\.csv$", full.names = TRUE)
  
  # Safety check: Ensure files were found
  if (length(mtx_file) == 0) stop("Error: No .mtx file found in ", data_dir)
  if (length(csv_file) == 0) stop("Error: No coordinate .csv file found in ", data_dir)
  
  message("Loading matrix: ", basename(mtx_file[1]))
  message("Loading coords: ", basename(csv_file[1]))
  
  # --- 2. Load and Prepare Data ---
  # Note: readMM creates a sparse matrix but doesn't include gene/cell names
  counts <- Matrix::readMM(mtx_file[1])
  bead_coords <- read.csv(csv_file[1], header = TRUE, check.names = FALSE)
  rownames(bead_coords) <- bead_coords[,1]
  bead_coords <- bead_coords[,-1]
  
  # Important: If your .mtx doesn't have names, you usually need genes.tsv/barcodes.tsv
  # We assume the matrix dimensions match the coordinate rows
  if (ncol(counts) != nrow(bead_coords)) {
    stop("Dimension mismatch: Matrix columns do not match coordinate rows.")
  }
  
  colnames(counts) <- rownames(bead_coords)
  
  # --- 3. Create Seurat Object ---
  seu <- CreateSeuratObject(counts = counts, assay = "Spatial")
  
  # --- 4. Format Spatial Metadata ---
  # To match your specific slide format: tissue, row, col, imagerow, imagecol
  spatial_coords <- data.frame(
    tissue = 1,                 
    row = bead_coords[[1]],     
    col = bead_coords[[2]],     
    imagerow = bead_coords[[1]],
    imagecol = bead_coords[[2]],
    row.names = rownames(bead_coords)
  )
  
  # --- 5. Build SlideSeq Image Slot ---
  # SlideSeq class is designed for non-gridded spatial data (beads)
  slide_image <- new(
    Class = 'SlideSeq',
    coordinates = spatial_coords,
    assay = "Spatial",
    key = "slice1_"
  )
  
  # Add to object
  seu[["slice1"]] <- slide_image
} else if (platform == "StereoSeq") {
  library(hdf5r)
  h5ad_files <- list.files(data_dir, pattern = "\\.h5ad$", full.names = TRUE)
  if (length(h5ad_files) == 0) stop("No .h5ad found for StereoSeq.")
  
  # Prioritize 'processed' if it exists, otherwise take first
  h5_path <- grep("_processed.h5ad$", h5ad_files, value = TRUE)[1]
  if (is.na(h5_path)) h5_path <- h5ad_files[1]
  
  message("Manually extracting data from H5AD: ", h5_path)
  file.h5 <- H5File$new(h5_path, mode = "r")
  
  # 1. DYNAMIC DATA PATH DETECTION
  counts_path <- NULL
  for (path in c("layers/raw_count", "X", "raw/X")) {
    if (file.h5$exists(path)) { counts_path <- path; break }
  }
  if (is.null(counts_path)) stop("Could not find counts in layers/raw_count, X, or raw/X.")
  
  # 2. DIMENSIONS & NAMES
  # Standard H5AD index paths
  obs_names <- if(file.h5$exists("obs/_index")) file.h5[["obs/_index"]][] else stop("Missing obs index")
  var_names <- if(file.h5$exists("var/_index")) file.h5[["var/_index"]][] else stop("Missing var index")
  
  n_obs <- length(obs_names)
  n_var <- length(var_names)
  message(sprintf(" - Found %d cells and %d genes", n_obs, n_var))
  
  # 3. SPARSE MATRIX RECONSTRUCTION
  if (is(file.h5[[counts_path]], "H5Group")) {
    message(" - Reconstruction from Sparse Group: ", counts_path)
    data <- file.h5[[paste0(counts_path, "/data")]][]
    indices <- file.h5[[paste0(counts_path, "/indices")]][]
    indptr <- file.h5[[paste0(counts_path, "/indptr")]][]
    
    # DETERMINE IF CSR OR CSC
    # If indptr length is n_obs + 1, it's CSR (Pointer = Cells, Indices = Genes)
    if (length(indptr) == (n_obs + 1)) {
      message("   - Matrix format: CSR (Cells as pointers)")
      # We want (Genes x Cells). By using p = indptr (Cells), 
      # indices + 1 must be rows (Genes).
      counts <- sparseMatrix(
        i = indices + 1,
        p = indptr,
        x = as.numeric(data),
        dims = c(n_var, n_obs),
        dimnames = list(var_names, obs_names),
        giveCsparse = TRUE
      )
    } else {
      message("   - Matrix format: CSC (Genes as pointers)")
      # pointer is Genes, so indices + 1 are Cells. This creates (Cells x Genes)
      counts_tmp <- sparseMatrix(
        i = indices + 1,
        p = indptr,
        x = as.numeric(data),
        dims = c(n_obs, n_var),
        dimnames = list(obs_names, var_names),
        giveCsparse = TRUE
      )
      counts <- Matrix::t(counts_tmp)
    }
  } else {
    message(" - Loading Dense Matrix...")
    counts <- t(file.h5[[counts_path]][,])
    rownames(counts) <- var_names
    colnames(counts) <- obs_names
  }
  
  # 4. COORDINATES
  message(" - Extracting spatial coordinates...")
  if (file.h5$exists("obs/x") && file.h5$exists("obs/y")) {
    x_coords <- file.h5[["obs/x"]][]
    y_coords <- file.h5[["obs/y"]][]
  } else if (file.h5$exists("obsm/spatial")) {
    # Scanpy obsm is stored as (n_obs, 2)
    spatial_mat <- file.h5[["obsm/spatial"]][,]
    # In HDF5r, matrices can be transposed depending on how they were written
    if (ncol(spatial_mat) == n_obs) {
      x_coords <- spatial_mat[1,]
      y_coords <- spatial_mat[2,]
    } else {
      x_coords <- spatial_mat[,1]
      y_coords <- spatial_mat[,2]
    }
  } else {
    stop("Coordinates not found in obs/x,y or obsm/spatial.")
  }
  
  # 5. ASSEMBLE
  seu <- CreateSeuratObject(counts = counts, assay = "Spatial")
  
  # Ensure any metadata from 'obs' is carried over
  obs_data <- as.data.frame(lapply(names(file.h5[["obs"]]), function(n) {
    if(is(file.h5[[paste0("obs/", n)]], "H5Dataset")) return(file.h5[[paste0("obs/", n)]][])
  }))
  # Filter out list elements that are NULL or complex
  obs_data <- obs_data[!sapply(obs_data, is.null)]
  
  spatial_coords <- data.frame(
    imagecol = x_coords,
    imagerow = y_coords,
    row = y_coords,
    col = x_coords,
    row.names = colnames(seu)
  )
  
  seu[["slice1"]] <- new(Class = "SlideSeq", assay = "Spatial", 
                         coordinates = spatial_coords, key = "slice1_")
  
  file.h5$close_all()
  message("Successfully assembled StereoSeq object.")
}

# Final standardization
if (!"Spatial" %in% names(seu@assays)) {
  hd_assays <- grep("^Spatial\\.", names(seu@assays), value = TRUE)
  if (length(hd_assays) > 0) {
    seu[["Spatial"]] <- seu[[hd_assays[1]]]
    DefaultAssay(seu) <- "Spatial"
    seu[[hd_assays[1]]] <- NULL
  }
} else {
  DefaultAssay(seu) <- "Spatial"
}

dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
saveRDS(seu, file = outfile)
message("Saved Seurat object to: ", outfile)
