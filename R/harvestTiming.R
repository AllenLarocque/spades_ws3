#' Apply one year's slice of a ws3 harvest plan to an inventory raster
#'
#' Used when `harvestTiming = "annual"`. ws3's spatial allocation (SDA) writes, in each
#' planning year, the inventory with the whole period's harvest already applied (`postPlan`)
#' and one planned-harvest raster per year of the period. This applies only the pixels planned
#' for the current year, taking their development type (band 1) and age (band 2) from
#' `postPlan`. A planned pixel whose age is now lower than its age at planning time plus the
#' years elapsed was reset since the plan (e.g. burned); it is skipped and counted.
#'
#' @param current  SpatRaster (3 bands: dev-type hash, age, block). Inventory before this year's harvest.
#' @param prePlan  SpatRaster. Inventory at the planning year, before the SDA.
#' @param postPlan SpatRaster. Inventory the SDA wrote at the planning year (whole period harvested).
#' @param harvest  SpatRaster. This year's planned harvest (`projected_harvest_<year>.tif`, 1 = harvest).
#' @param elapsed  Integer. Years since the planning year.
#' @return List: `inventory` (SpatRaster, in memory), `applied_mask` (SpatRaster, 1 = harvested,
#'   NA otherwise), `planned`, `applied`, `skippedReset` (pixel counts).
applyAnnualHarvest <- function(current, prePlan, postPlan, harvest, elapsed) {
  cur  <- terra::values(current)
  pre  <- terra::values(prePlan)
  post <- terra::values(postPlan)
  h    <- terra::values(harvest)[, 1]
  planned <- !is.na(h) & h == 1
  reset   <- planned & !is.na(cur[, 2]) & !is.na(pre[, 2]) & cur[, 2] < pre[, 2] + elapsed
  apply   <- planned & !reset
  cur[apply, 1:2] <- post[apply, 1:2]
  out <- terra::rast(current)
  terra::values(out) <- cur
  names(out) <- names(current)
  appliedMask <- terra::rast(current, nlyrs = 1)
  terra::values(appliedMask) <- ifelse(apply, 1L, NA_integer_)
  list(inventory = out, applied_mask = appliedMask,
       planned = sum(planned), applied = sum(apply), skippedReset = sum(reset))
}
