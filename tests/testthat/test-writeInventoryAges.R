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

test_that("applyHarvest copies age resets into this year's inventory before ws3 reads it", {
  # Regression: it wrote next year's file (offset = 1), so burned-age resets never reached ws3
  src <- readLines(test_path("..", "..", "spades_ws3.R"))
  body <- paste(src[grep("^applyHarvest <- function", src):grep("py\\$simulate_harvest", src)[1]],
                collapse = "\n")
  expect_match(body, "syncAgeResets\\(sim\\$landscape\\$age, inventoryFilesForYear\\(sim, year\\)\\)")
})

## Two overlapping management units on one grid: unit A covers columns 1-3, unit B columns 2-4
twoUnits <- function(d) {
  base <- terra::rast(nrows = 1, ncols = 4, xmin = 0, xmax = 4, ymin = 0, ymax = 1)
  mk <- function(cols, ages, f) {
    r <- terra::rast(base, nlyrs = 3)
    a <- rep(NA_real_, 4); a[cols] <- ages
    terra::values(r) <- cbind(ifelse(is.na(a), NA, 1), a, ifelse(is.na(a), NA, 1))
    terra::writeRaster(terra::crop(r, terra::ext(min(cols) - 1, max(cols), 0, 1)), f,
                       filetype = "GTiff", datatype = "INT4S", overwrite = TRUE)
    f
  }
  list(base = base,
       A = mk(1:3, c(50, 80, 1), file.path(d, "A.tif")),   # unit A harvested column 3 (age 1)
       B = mk(2:4, c(60, 90, 70), file.path(d, "B.tif")))  # overlap: A says 80/1, B says 60/90
}

test_that("syncAgeResets keeps each unit's own ages and copies only resets to 0", {
  d <- withr::local_tempdir()
  u <- twoUnits(d)
  landscape <- terra::rast(u$base); terra::values(landscape) <- c(50, 70, 0, 70)  # merged ages; column 3 burned
  n <- syncAgeResets(landscape, c(u$A, u$B))
  expect_equal(n, c(1L, 1L))                                         # the shared burned cell, in both units
  expect_equal(terra::values(terra::rast(u$A))[, 2], c(50, 80, 0))   # column 2 keeps A's 80, not the merged 70
  expect_equal(terra::values(terra::rast(u$B))[, 2], c(60, 0, 70))   # column 2 keeps B's 60
})

test_that("syncAgeResets writes nothing when there are no resets (standalone runs)", {
  d <- withr::local_tempdir()
  u <- twoUnits(d)
  before <- tools::md5sum(c(u$A, u$B))
  landscape <- terra::rast(u$base); terra::values(landscape) <- c(50, 70, 1, 70)  # merged view, no zeros
  expect_equal(syncAgeResets(landscape, c(u$A, u$B)), c(0L, 0L))
  expect_identical(tools::md5sum(c(u$A, u$B)), before)
})

test_that("syncAgeResets does not touch inventory cells that are already 0 or NA", {
  d <- withr::local_tempdir()
  f <- makeInventory(file.path(d, "inv.tif"), c(0, 60, 70, 80, NA, 90))
  age <- terra::rast(terra::rast(f), nlyrs = 1)
  terra::values(age) <- c(0, 60, 0, 80, 0, 90)
  expect_equal(syncAgeResets(age, f), 1L)
  expect_equal(terra::values(terra::rast(f))[, 2], c(0, 60, 0, 80, NA, 90))
})
