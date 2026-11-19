#' Write stand ages into ws3 inventory rasters
#'
#' Replaces the age band (band 2) of each source inventory raster with `age`, cropped and
#' masked to that raster, and writes the result to the matching destination file. Other
#' bands are copied unchanged. A destination may be its own source (the raster is read into
#' memory first).
#'
#' @param age      SpatRaster of stand ages covering all management units (`sim$landscape$age`).
#' @param srcFiles Character. Inventory GeoTIFFs to read (one per management unit).
#' @param dstFiles Character. GeoTIFFs to write, in the same order as `srcFiles`.
#' @return `dstFiles`, invisibly.
writeInventoryAges <- function(age, srcFiles, dstFiles) {
  stopifnot(length(srcFiles) == length(dstFiles))
  for (i in seq_along(srcFiles)) {
    rs <- terra::toMemory(terra::rast(srcFiles[i]))  # in memory, so dstFiles[i] may equal srcFiles[i]
    newAge <- terra::mask(terra::crop(age, rs[[2]]), rs[[2]])
    newAge[is.nan(newAge)] <- NA
    rs[[2]] <- newAge
    terra::writeRaster(rs, dstFiles[i], filetype = "GTiff", overwrite = TRUE, datatype = "INT4S")
  }
  invisible(dstFiles)
}
