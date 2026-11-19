# transform VisiumV2-style image to VisiumV1-style
# input:
#   path_to_visium : 路径（示例 "/.../Mouse_embryo_Visium/"），函数会在此路径下查找 spatial/tissue_positions.csv
#   seu             : 已存在的 Seurat 对象（将被拷贝 image 到返回对象）
#   image_name      : image slice 名称（默认 "slice1"）
# output:
#   返回修改后的 Seurat 对象（内部将用 VisiumV1 替换 seu[[image_name]]）
convert_visium_v2_to_v1 <- function(path_to_visium, seu, image_name = "slice1") {
  library(Seurat)
  
  path_to_coordinates <- file.path(path_to_visium, "spatial", "tissue_positions.csv")
  object <- Seurat::Load10X_Spatial(path_to_visium, slice = image_name)
  visium_v2 <- object[[image_name]]
  coordinates <- Seurat::Read10X_Coordinates(path_to_coordinates, filter.matrix = TRUE)
  
  visium_v1 <- new(
    Class = "VisiumV1",
    assay = visium_v2@assay,
    key = visium_v2@key,
    coordinates = coordinates,
    scale.factors = visium_v2@scale.factors,
    image = visium_v2@image
  )
  
  visium_v1@spot.radius <- Seurat::Radius(visium_v1)
  
  # align counts & coordinates order
  visium_v1 <- visium_v1[Cells(seu)]
  seu[[image_name]] <- visium_v1
  
  return(seu)
}

seu_fixed <- convert_visium_v2_to_v1("/vast/projects/SpatialTranscriptomicsBenchmarking/Data/Mouse_embryo_Visium/",
                                     seu,
                                     image_name = "slice1")
