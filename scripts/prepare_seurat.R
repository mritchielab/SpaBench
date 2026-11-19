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
  make_option(c("--platform"), type="character", default="Visium", help="Platform: Visium, VisiumHD, Slideseq, StereoSeq"),
  make_option(c("--bin.size"), type="numeric", default=8, help="Bin size for VisiumHD: 8 or 16"),
  # ---- Step 0 QC ----
  make_option(c("--remove_mt"), type="logical", default=TRUE,
              help="QC: remove mitochondrial genes (all platforms) [default TRUE]"),
  make_option(c("--mt_pattern"), type="character", default="^mt-",
              help="QC: regular expression for mitochondrial gene names, matched case-insensitively [default '^mt-']"),
  make_option(c("--in_tissue_only"), type="logical", default=TRUE,
              help="QC: Visium / Visium HD - keep only spots/bins with in_tissue = 1 [default TRUE]")
)

opt <- parse_args(OptionParser(option_list=option_list))

data_dir <- opt$data_dir
filename <- opt$filename
outfile  <- opt$out
platform <- opt$platform
# Accept any capitalisation of the platform name and map it to the spelling used below
platform <- switch(tolower(platform), visium = "Visium", visiumhd = "VisiumHD",
                   stereoseq = "StereoSeq", slideseq = "SlideSeq", platform)
bin_size <- opt$bin.size

if (!dir.exists(data_dir)) stop("data_dir not found: ", data_dir)

# ==============================================================================
# Step 0 QC (same rules in scripts/prepare_anndata.py)
#   * mitochondrial genes are removed for every platform
#   * Visium / Visium HD: spots / bins outside the captured tissue (in_tissue != 1
#     in tissue_positions) are removed
# ==============================================================================
remove_mt      <- isTRUE(opt$remove_mt)
in_tissue_only <- isTRUE(opt$in_tissue_only)
qc_remove_mt <- function(counts) {
  if (!remove_mt) { message("QC: mitochondrial genes kept (--remove_mt FALSE)"); return(counts) }
  mt <- grepl(opt$mt_pattern, rownames(counts), ignore.case = TRUE)
  message(sprintf("QC: removing %d mitochondrial genes (pattern '%s', case-insensitive)%s",
                  sum(mt), opt$mt_pattern,
                  if (any(mt)) paste0(": ", paste(head(rownames(counts)[mt], 15), collapse = ", ")) else ""))
  if (!any(mt)) warning("QC: no gene matched the mitochondrial pattern '", opt$mt_pattern, "'")
  counts[!mt, , drop = FALSE]
}

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
  if (is.list(counts)) counts <- counts[["Gene Expression"]]   # multi-modal .h5
  counts <- qc_remove_mt(counts)
  # Ensure the assay is named "Spatial" from the start
  seu <- CreateSeuratObject(counts = counts, assay = "Spatial")
  tissue_filtered <- FALSE
  
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
        # Space Ranger >= 2.0 writes tissue_positions.csv WITH a header row; older
        # tissue_positions_list.csv has none. Detect it so the header is not read as data.
        has_header <- grepl("barcode", readLines(coord_file, n = 1L), ignore.case = TRUE)
        positions <- read.csv(coord_file, header = has_header, row.names = 1)
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
      
      # positions files without the "-1" barcode suffix used in the matrix
      if (!any(colnames(seu) %in% rownames(positions)) &&
          any(colnames(seu) %in% paste0(rownames(positions), "-1"))) {
        rownames(positions) <- paste0(rownames(positions), "-1")
      }

      # QC: keep only spots / bins inside the captured tissue
      if (in_tissue_only) {
        tissue_bc <- rownames(positions)[as.integer(positions$tissue) == 1]
        keep <- colnames(seu)[colnames(seu) %in% tissue_bc]
        if (length(keep) == 0) stop("QC: no matrix barcode has in_tissue = 1 in ", coord_file)
        message(sprintf("QC: keeping %d of %d spots/bins inside the tissue (in_tissue = 1 in %s)",
                        length(keep), ncol(seu), basename(coord_file)))
        if (length(keep) < ncol(seu)) seu <- subset(seu, cells = keep)
        tissue_filtered <- TRUE
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
  if (in_tissue_only && !tissue_filtered)
    stop("QC: no tissue_positions file found for ", platform, " in ", data_dir,
         " - cannot remove spots outside the tissue (set qc: in_tissue_only: false to skip)")
} else if (platform %in% c("StereoSeq", "SlideSeq")) {
  # ----------------------------------------------------------------------------
  # Slide-seq / Stereo-seq. File discovery, count selection and bead/bin
  # filtering mirror scripts/prepare_anndata.py so that 01_raw_seurat.rds and
  # the two 01_*anndata.h5ad files hold the same barcodes and genes.
  # Files written by the workflow itself (01_*) are never used as input.
  # ----------------------------------------------------------------------------
  library(hdf5r)

  list_inputs <- function(d, pattern) {
    f <- list.files(d, pattern = pattern, ignore.case = TRUE)
    f <- f[!startsWith(f, "01_") & !dir.exists(file.path(d, f))]
    sort(f, method = "radix")            # C-locale order, same as Python sorted()
  }
  is_integer_valued <- function(x) {
    length(x) == 0 || (all(is.finite(x)) && all(x >= 0) && all(x == round(x)))
  }
  # h5$exists() errors when an intermediate group (e.g. "raw") is missing
  h5_has <- function(h5, path) {
    parts <- strsplit(path, "/", fixed = TRUE)[[1]]
    for (k in seq_along(parts)) {
      if (!h5$exists(paste(parts[1:k], collapse = "/"))) return(FALSE)
    }
    TRUE
  }
  h5_index <- function(h5, grp) {
    col <- if ("_index" %in% h5attr_names(h5[[grp]])) h5attr(h5[[grp]], "_index") else "_index"
    as.character(h5[[paste0(grp, "/", col)]][])
  }
  # AnnData matrix element -> genes x cells dgCMatrix (h5ad stores cells x genes)
  read_h5_matrix <- function(h5, key, n_obs, n_var) {
    obj <- h5[[key]]
    if (inherits(obj, "H5Group")) {
      enc <- h5attr(obj, "encoding-type")
      x   <- as.numeric(obj[["data"]][])
      i   <- as.integer(obj[["indices"]][]) + 1L
      p   <- as.integer(obj[["indptr"]][])
      if (enc == "csr_matrix") {         # pointers = cells -> columns
        sparseMatrix(i = i, p = p, x = x, dims = c(n_var, n_obs))
      } else {                           # csc: pointers = genes
        Matrix::t(sparseMatrix(i = i, p = p, x = x, dims = c(n_obs, n_var)))
      }
    } else {
      m <- obj[, ]                       # hdf5r returns the transposed (genes x cells) view
      if (!(nrow(m) == n_var && ncol(m) == n_obs)) m <- t(m)
      as(Matrix(m, sparse = TRUE), "CsparseMatrix")
    }
  }
  count_keys <- c("layers/raw_count", "layers/count", "layers/counts", "raw/X", "X")
  # Raw integer counts (genes x cells), names and x/y coordinates from an .h5ad
  read_h5ad_counts <- function(path, need_coords = TRUE) {
    h5 <- H5File$new(path, mode = "r")
    on.exit(h5$close_all())
    obs_names <- h5_index(h5, "obs")
    counts <- NULL
    for (k in count_keys) {
      if (!h5_has(h5, k)) next
      vgrp <- if (startsWith(k, "raw/")) "raw/var" else "var"
      var_names <- h5_index(h5, vgrp)
      m <- read_h5_matrix(h5, k, length(obs_names), length(var_names))
      if (is_integer_valued(m@x)) { counts <- m; key <- k; break }
      message("  ", k, ": not integer-valued (normalised?) - skipped")
    }
    if (is.null(counts))
      stop(path, ": no integer-valued count matrix in ", paste(count_keys, collapse = ", "),
           "; raw counts are required (only normalised values found).")
    dimnames(counts) <- list(var_names, obs_names)
    x <- y <- NULL
    if (need_coords) {
      if (h5_has(h5, "obs/x") && h5_has(h5, "obs/y")) {
        x <- as.numeric(h5[["obs/x"]][]); y <- as.numeric(h5[["obs/y"]][])
      } else if (h5_has(h5, "obsm/spatial")) {
        s <- h5[["obsm/spatial"]][, ]
        if (ncol(s) == length(obs_names) && nrow(s) != length(obs_names)) s <- t(s)
        x <- as.numeric(s[, 1]); y <- as.numeric(s[, 2])
      } else stop(path, ": coordinates not found in obs/x,y or obsm/spatial.")
    }
    message(sprintf("  counts from '%s': %d cells x %d genes, total = %.0f",
                    key, ncol(counts), nrow(counts), sum(counts@x)))
    list(counts = counts, x = x, y = y, key = key)
  }
  read_mtx <- function(path) {
    con <- if (grepl("\\.gz$", path)) gzfile(path, "r") else file(path, "r")
    hdr <- readLines(con, 1); skip <- 1L
    repeat { l <- readLines(con, 1); skip <- skip + 1L; if (!startsWith(l, "%")) break }
    close(con)
    if (!grepl("coordinate", hdr) || !grepl("general", hdr)) {
      return(as(Matrix::readMM(path), "CsparseMatrix"))
    }
    dims <- as.numeric(strsplit(trimws(l), "\\s+")[[1]])
    dt <- data.table::fread(path, skip = skip, header = FALSE, sep = " ",
                            col.names = c("i", "j", "x"))
    sparseMatrix(i = dt$i, j = dt$j, x = as.numeric(dt$x), dims = dims[1:2])
  }
  # first column (and header) of a dense genes x beads digital_expression.txt
  read_dge_names <- function(path) {
    con <- if (grepl("\\.gz$", path)) gzfile(path, "r") else file(path, "r")
    on.exit(close(con))
    bc <- strsplit(readLines(con, 1), "\t", fixed = TRUE)[[1]][-1]
    genes <- character(0)
    repeat {
      l <- readLines(con, n = 500)
      if (length(l) == 0) break
      genes <- c(genes, sub("\t.*$", "", l))
    }
    list(barcodes = bc, genes = genes)
  }
  read_first_col <- function(path) {
    as.character(data.table::fread(path, header = FALSE, select = 1)[[1]])
  }

  if (platform == "StereoSeq") {
    h5ad_files <- list_inputs(data_dir, "\\.h5ad$")
    if (length(h5ad_files) == 0) stop("No input .h5ad (other than 01_*) found for StereoSeq in ", data_dir)
    # Prefer *_processed.h5ad, then the --filename passed by the Snakefile, else the first
    proc <- h5ad_files[grepl("_processed\\.h5ad$", h5ad_files, ignore.case = TRUE)]
    h5_name <- if (length(proc) > 0) proc[1] else if (filename %in% h5ad_files) filename else h5ad_files[1]
    h5_path <- file.path(data_dir, h5_name)
    message("Reading Stereo-seq counts from H5AD: ", h5_path)
    r <- read_h5ad_counts(h5_path)
    counts <- r$counts; x_coords <- r$x; y_coords <- r$y

  } else {  # SlideSeq
    # Built ONLY from the raw bead files: *expression_matrix.mtx + *locations*.csv
    # (barcodes / gene names from a barcodes / genes file or digital_expression.txt).
    # Any .h5ad in the folder is ignored. All beads present in both the matrix
    # and the locations file are kept.
    mtx  <- list_inputs(data_dir, "expression_matrix\\.mtx(\\.gz)?$")
    locs <- list_inputs(data_dir, "locations.*\\.csv(\\.gz)?$")
    if (length(mtx) == 0 || length(locs) == 0)
      stop("Slide-seq input needs *expression_matrix.mtx and *locations*.csv in ", data_dir)
    mtx_path <- file.path(data_dir, mtx[1]); loc_path <- file.path(data_dir, locs[1])
    message("  matrix:    ", mtx_path)
    message("  locations: ", loc_path)
    M <- read_mtx(mtx_path)

    loc <- read.csv(loc_path, check.names = FALSE, stringsAsFactors = FALSE)
    loc_bc <- trimws(as.character(loc[[1]]))
    keep_loc <- !duplicated(loc_bc)
    loc <- loc[keep_loc, , drop = FALSE]; rownames(loc) <- loc_bc[keep_loc]
    lc <- tolower(colnames(loc))
    xcol <- if ("xcoord" %in% lc) which(lc == "xcoord")[1] else if ("x" %in% lc) which(lc == "x")[1] else 2
    ycol <- if ("ycoord" %in% lc) which(lc == "ycoord")[1] else if ("y" %in% lc) which(lc == "y")[1] else 3

    # --- bead barcodes (matrix order) ---
    bc_files <- list_inputs(data_dir, "barcodes.*\\.(tsv|txt|csv)(\\.gz)?$")
    dge      <- list_inputs(data_dir, "digital_expression\\.txt(\\.gz)?$")
    dge_names <- NULL
    if (length(bc_files) > 0) {
      barcodes <- trimws(read_first_col(file.path(data_dir, bc_files[1]))); bc_src <- bc_files[1]
    } else if (length(dge) > 0) {
      dge_names <- read_dge_names(file.path(data_dir, dge[1]))
      barcodes <- dge_names$barcodes; bc_src <- paste0(dge[1], " (header)")
    } else {
      barcodes <- rownames(loc); bc_src <- paste0(locs[1], " (row order)")
    }
    if (nrow(M) == length(barcodes) && ncol(M) != length(barcodes)) {
      M <- Matrix::t(M)                  # beads x genes -> genes x beads
    } else if (ncol(M) != length(barcodes)) {
      stop(sprintf("matrix %d x %d matches neither dimension of %d barcodes from %s",
                   nrow(M), ncol(M), length(barcodes), bc_src))
    }
    message("  barcodes from ", bc_src, ": ", length(barcodes))
    n_genes <- nrow(M)

    # --- gene names ---
    genes <- NULL; g_src <- NULL
    g_files <- list_inputs(data_dir, "(genes|features).*\\.(tsv|txt|csv)(\\.gz)?$")
    if (length(g_files) > 0) {
      g <- read_first_col(file.path(data_dir, g_files[1]))
      if (length(g) == n_genes) { genes <- g; g_src <- g_files[1] }
    }
    if (is.null(genes) && length(dge) > 0) {
      if (is.null(dge_names)) dge_names <- read_dge_names(file.path(data_dir, dge[1]))
      if (length(dge_names$genes) == n_genes) {
        genes <- dge_names$genes; g_src <- paste0(dge[1], " (first column)")
      }
    }
    if (is.null(genes)) {
      warning("No gene names found for the Slide-seq matrix (no genes/features file ",
              "or digital_expression.txt) - using placeholders Feature1..N")
      genes <- paste0("Feature", seq_len(n_genes)); g_src <- "PLACEHOLDER"
    }
    message("  genes from ", g_src, ": ", n_genes)

    # --- keep beads present in both the matrix and the locations file ---
    keep <- which(barcodes %in% rownames(loc))
    message(sprintf("  beads: matrix %d, locations %d, kept (both) %d",
                    length(barcodes), nrow(loc), length(keep)))
    if (length(keep) == 0) stop("no bead barcodes shared between matrix and locations file")
    counts <- M[, keep, drop = FALSE]
    dimnames(counts) <- list(genes, barcodes[keep])
    loc <- loc[barcodes[keep], , drop = FALSE]
    x_coords <- as.numeric(loc[[xcol]]); y_coords <- as.numeric(loc[[ycol]])
  }

  # Gene names: Seurat turns "_" into "-"; do it explicitly and de-duplicate the
  # same way as prepare_anndata.py so both objects carry identical names.
  rownames(counts) <- make.unique(gsub("_", "-", rownames(counts), fixed = TRUE), sep = "-")
  counts <- as(counts, "CsparseMatrix")
  counts <- qc_remove_mt(counts)   # Step 0 QC (Slide-seq / Stereo-seq have no off-tissue flag)

  seu <- CreateSeuratObject(counts = counts, assay = "Spatial")

  # Image coordinates: x/y first (Seurat's SlideSeq accessors use the first two
  # columns), plus row/col and imagerow/imagecol used by the downstream scripts.
  spatial_coords <- data.frame(
    x = x_coords, y = y_coords,
    imagerow = y_coords, imagecol = x_coords,
    row = y_coords, col = x_coords,
    row.names = colnames(seu)
  )
  seu[["slice1"]] <- new(Class = "SlideSeq", assay = "Spatial",
                         coordinates = spatial_coords, key = "slice1_")
  message("Successfully assembled ", platform, " object: ",
          nrow(seu), " genes x ", ncol(seu), " spots.")
} else {
  stop("Unknown platform '", opt$platform, "' (expected Visium, VisiumHD, StereoSeq or SlideSeq)")
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