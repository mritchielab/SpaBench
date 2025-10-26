#!/usr/bin/env Rscript
# scripts/prepare_seurat.R
# Usage:
# Rscript scripts/prepare_seurat.R --data_dir ./path/ --filename filtered_feature_bc_matrix.h5 --out data/SAMPLE/01_raw_seurat.rds --platform Visium

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
})

option_list <- list(
  make_option(c("--data_dir"), type="character", help="Directory containing input files"),
  make_option(c("--filename"), type="character", default="filtered_feature_bc_matrix.h5", help="HDF5 filename (for Visium)"),
  make_option(c("--out"), type="character", help="Output .rds path"),
  make_option(c("--platform"), type="character", default="Visium", help="Platform: Visium, VisiumHD, SlideSeq, StereoSeq")
)

opt <- parse_args(OptionParser(option_list=option_list))

data_dir <- opt$data_dir
filename <- opt$filename
outfile  <- opt$out
platform <- opt$platform

if (!dir.exists(data_dir)) {
  stop("data_dir not found: ", data_dir)
}

# Only check h5 file for Visium/VisiumHD to preserve original behavior
if (platform %in% c("Visium", "VisiumHD")) {
  h5path <- file.path(data_dir, filename)
  if (!file.exists(h5path)) stop("H5 file not found: ", h5path)
}

message(Sys.time(), " - Creating Seurat object for platform: ", platform)

# --------- minimal platform-specific creation logic ----------
if (platform %in% c("Visium", "VisiumHD")) {
  # keep original creation logic
  seu <- tryCatch({
    Load10X_Spatial(data.dir = data_dir, filename = filename, assay = "Spatial", slice = "slice1")
  }, error = function(e) stop("Load10X_Spatial failed: ", e$message))
  
} else if (platform == "SlideSeq") {
  # minimal SlideSeq handling: try Read10X_h5() then Read10X(); fail otherwise
  h5file <- file.path(data_dir, filename)
  counts <- NULL
  if (file.exists(h5file)) {
    counts <- tryCatch(Read10X_h5(h5file), error = function(e) { message("Read10X_h5 failed: ", e$message); NULL })
    if (is.list(counts)) counts <- counts[[1]]
  }
  if (is.null(counts)) {
    counts <- tryCatch(Read10X(data.dir = data_dir), error = function(e) { stop("Read10X(data.dir) failed: ", e$message) })
    if (is.list(counts)) counts <- counts[[1]]
  }
  # Create Seurat object from counts (assay = Spatial) — mirrors original downstream assumptions
  seu <- CreateSeuratObject(counts = counts, assay = "Spatial", project = basename(data_dir))
  
} else if (platform == "StereoSeq") {
  # prefer .h5seurat, else .h5ad convert -> h5seurat
  h5seurat_files <- list.files(data_dir, pattern = "\\.h5seurat$", full.names = TRUE)
  if (length(h5seurat_files) >= 1) {
    h5s <- h5seurat_files[1]
    seu <- tryCatch({
      SeuratDisk::LoadH5Seurat(h5s, assays = "Spatial", images = FALSE)
    }, error = function(e) stop("LoadH5Seurat failed: ", e$message))
  } else {
    h5ad_files <- list.files(data_dir, pattern = "\\.h5ad$", full.names = TRUE)
    if (length(h5ad_files) >= 1) {
      h5ad <- h5ad_files[1]
      SeuratDisk::Convert(h5ad, dest = "h5seurat", overwrite=TRUE, verbose=TRUE, assay="Spatial")
      h5s_out <- file.path(dirname(h5ad), paste0(tools::file_path_sans_ext(basename(h5ad)), ".h5seurat"))
      seu <- tryCatch({
        SeuratDisk::LoadH5Seurat(h5s_out, assays = "Spatial", images = FALSE)
      }, error = function(e) stop("LoadH5Seurat (after convert) failed: ", e$message))
    } else {
      stop("StereoSeq selected but no .h5seurat or .h5ad found in data_dir.")
    }
  }
  
  ## --- 统一在加载后解析 barcodes -> coords 并插入 image（或回退写 meta.data） ---
  barcodes <- rownames(seu@meta.data)
  nums_list <- regmatches(barcodes, gregexpr("\\d+", barcodes))
  coords_mat <- do.call(rbind, lapply(nums_list, function(v) {
    v <- as.numeric(v)
    if (length(v) >= 2) return(v[1:2]) else return(c(NA, NA))
  }))
  coords <- as.data.frame(coords_mat, stringsAsFactors = FALSE)
  colnames(coords) <- c("imagecol", "imagerow")
  coords$row <- coords$imagerow
  coords$col <- coords$imagecol
  rownames(coords) <- barcodes
  
  # 尝试创建 SlideSeq image 并放入 seu@images；若失败则将 coords 写入 meta.data 回退
  safe_create_image <- function(seu, coords_df, image_name = "slice1") {
    ok <- FALSE
    result_seu <- seu
    try({
      img_obj <- new(Class = "SlideSeq", assay = "Spatial", coordinates = coords_df)
      result_seu@images[[image_name]] <- img_obj
      # try to set key if slot exists
      if (isS4(result_seu@images[[image_name]]) && "key" %in% slotNames(class(result_seu@images[[image_name]]))) {
        slot(result_seu@images[[image_name]], "key") <- paste0(image_name, "_")
      }
      ok <- TRUE
    }, silent = TRUE)
    if (!ok) {
      # fallback: merge coords into meta.data
      result_seu@meta.data <- cbind(result_seu@meta.data, coords_df[, c("col","row")])
      message("Note: SlideSeq image class creation failed; coords were added to seu@meta.data as 'col' and 'row'.")
    }
    return(result_seu)
  }
  
  seu <- safe_create_image(seu, coords, image_name = "slice1")
  if (!("nCount_Spatial" %in% colnames(seu@meta.data))) {
    if ("total_counts" %in% colnames(seu@meta.data)) {
      seu@meta.data$nCount_Spatial <- seu@meta.data$total_counts
    }
  }
  if (!("nFeature_Spatial" %in% colnames(seu@meta.data))) {
    if ("n_genes_by_counts" %in% colnames(seu@meta.data)) {
      seu@meta.data$nFeature_Spatial <- seu@meta.data$n_genes_by_counts
    }
  }
  
} else {
  stop("Unsupported platform: ", platform, ". Supported: Visium, VisiumHD, SlideSeq, StereoSeq")
}

# --------- keep original compatibility & save logic unchanged ----------
# Compatibility: VisiumHD often stores assay as Spatial.008um; rename/copy to 'Spatial'
if (! "Spatial" %in% names(seu@assays)) {
  if ("Spatial.008um" %in% names(seu@assays)) {
    message("Found 'Spatial.008um' assay -> copying to 'Spatial' and removing original")
    seu[["Spatial"]] <- seu[["Spatial.008um"]]
    DefaultAssay(seu) <- "Spatial"
    seu[["Spatial.008um"]] <- NULL
  } else {
    message("Warning: no 'Spatial' assay found after loading. assays: ", paste(names(seu@assays), collapse=", "))
  }
}

# Optionally: if platform is VisiumHD, we might want to remove very large image to save R memory.
# (Left commented — enable if necessary)
# if (platform == "VisiumHD") {
#   imgs <- names(seu@images)
#   if (length(imgs)>0) {
#     # remove raw image array to save memory: keep scale factors and coordinates
#     if (!is.null(seu@images[[imgs[1]]]@image)) {
#       seu@images[[imgs[1]]]@image <- NULL
#       message("Removed embedded image from Seurat object to save memory (VisiumHD).")
#     }
#   }
# }

dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
saveRDS(seu, file = outfile)
message("Saved Seurat object to: ", outfile)
