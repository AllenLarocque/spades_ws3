source(test_path("..", "..", "R", "inventoryAges.R"))

makeInventory <- function(path, ages) {
  r <- terra::rast(nrows = 2, ncols = 3, nlyrs = 3, vals = NA_real_)
  terra::values(r) <- cbind(c(11, 11, 12, 12, NA, 13),   # band 1: development-type hash
                            ages,                         # band 2: age
                            c(1, 1, 2, 2, NA, 3))         # band 3: block id
  terra::writeRaster(r, path, filetype = "GTiff", datatype = "INT4S", overwrite = TRUE)
  path
}

test_that("writeInventoryAges overwrites the age band of the file it reads (offset 0)", {
  d <- withr::local_tempdir()
  f <- makeInventory(file.path(d, "inventory_2030.tif"), c(50, 60, 70, 80, NA, 90))
  before <- terra::values(terra::rast(f))
  age <- terra::rast(terra::rast(f), nlyrs = 1)
  terra::values(age) <- c(51, 0, 71, 0, 5, 91)  # two burned pixels reset to 0; a value outside the inventory
  writeInventoryAges(age, f, f)
  after <- terra::values(terra::rast(f))
  expect_equal(after[, 2], c(51, 0, 71, 0, NA, 91))   # masked to the inventory: NA stays NA
  expect_equal(after[, c(1, 3)], before[, c(1, 3)])   # other bands unchanged
})

test_that("writeInventoryAges writes to a different destination without touching the source", {
  d <- withr::local_tempdir()
  src <- makeInventory(file.path(d, "inventory_2030.tif"), c(50, 60, 70, 80, NA, 90))
  dst <- file.path(d, "inventory_2031.tif")
  age <- terra::rast(terra::rast(src), nlyrs = 1)
  terra::values(age) <- c(51, 61, 71, 81, NA, 91)
  writeInventoryAges(age, src, dst)
  expect_equal(terra::values(terra::rast(dst))[, 2], c(51, 61, 71, 81, NA, 91))
  expect_equal(terra::values(terra::rast(src))[, 2], c(50, 60, 70, 80, NA, 90))
})

test_that("applyHarvest syncs landscape ages into this year's inventory before ws3 reads it", {
  # Regression: it wrote next year's file (offset = 1), so burned-age resets never reached ws3
  src <- readLines(test_path("..", "..", "spades_ws3.R"))
  body <- paste(src[grep("^applyHarvest <- function", src):grep("py\\$simulate_harvest", src)[1]],
                collapse = "\n")
  expect_match(body, "updateAges\\(sim, offset = 0\\)")
})
