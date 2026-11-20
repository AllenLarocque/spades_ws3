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

#' Copy age resets from the landscape into ws3 inventory rasters
#'
#' Other modules signal a stand-replacing disturbance (e.g. fire) by setting
#' `sim$landscape$age` to 0. This writes age 0 into each management unit's inventory where the
#' landscape age is 0 and the inventory age is above 0, and leaves every other age as the unit
#' has it (units may overlap, and the merged landscape does not preserve each unit's ages).
#' A reset on a cell shared by several units resets it in all of them. A file is only rewritten
#' when it has resets, so runs without such modules leave the inventory untouched.
#'
#' @param age      SpatRaster of stand ages covering all management units (`sim$landscape$age`).
#' @param invFiles Character. Inventory GeoTIFFs to update in place (one per management unit).
#' @return Integer vector: number of pixels reset in each file.
syncAgeResets <- function(age, invFiles) {
  vapply(invFiles, function(f) {
    rs <- terra::toMemory(terra::rast(f))
    landAge <- terra::values(terra::crop(age, rs[[2]]))[, 1]
    invAge  <- terra::values(rs[[2]])[, 1]
    if (length(landAge) != length(invAge))
      stop("syncAgeResets(): landscape and ", basename(f), " are not on the same grid")
    reset <- which(!is.na(landAge) & landAge == 0 & !is.na(invAge) & invAge > 0)
    if (length(reset)) {
      v <- terra::values(rs)
      v[reset, 2] <- 0
      terra::values(rs) <- v
      terra::writeRaster(rs, f, filetype = "GTiff", overwrite = TRUE, datatype = "INT4S")
    }
    length(reset)
  }, integer(1), USE.NAMES = FALSE)
}
