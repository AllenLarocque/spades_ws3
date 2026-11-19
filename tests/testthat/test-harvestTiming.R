source(test_path("..", "..", "R", "harvestTiming.R"))

inv <- function(dt, age, block = c(1, 1, 2, 2)) {
  r <- terra::rast(nrows = 2, ncols = 2, nlyrs = 3)
  terra::values(r) <- cbind(dt, age, block)
  r
}
one <- function(v) { r <- terra::rast(nrows = 2, ncols = 2); terra::values(r) <- v; r }

prePlan  <- inv(c(11, 11, 12, 12), c(80, 90, 100, 110))
postPlan <- inv(c(21, 21, 22, 12), c(0, 0, 0, 110))   # SDA harvested pixels 1-3 over the period

test_that("only this year's planned pixels are harvested, with type and age from the SDA", {
  current <- inv(c(11, 11, 12, 12), c(83, 93, 103, 113))   # 3 years after the plan
  res <- applyAnnualHarvest(current, prePlan, postPlan, one(c(NA, 1, NA, NA)), elapsed = 3)
  v <- terra::values(res$inventory)
  expect_equal(v[, 1], c(11, 21, 12, 12))
  expect_equal(v[, 2], c(83, 0, 103, 113))
  expect_equal(v[, 3], c(1, 1, 2, 2))
  expect_equal(c(res$planned, res$applied, res$skippedReset), c(1, 1, 0))
})

test_that("a planned pixel reset since the plan (burned) is skipped and counted", {
  current <- inv(c(11, 11, 12, 12), c(83, 2, 103, 113))     # pixel 2 burned a year ago
  res <- applyAnnualHarvest(current, prePlan, postPlan, one(c(1, 1, NA, NA)), elapsed = 3)
  v <- terra::values(res$inventory)
  expect_equal(v[, 2], c(0, 2, 103, 113))
  expect_equal(v[, 1], c(21, 11, 12, 12))
  expect_equal(c(res$planned, res$applied, res$skippedReset), c(2, 1, 1))
})

test_that("in the planning year the slice applies to the pre-plan inventory", {
  res <- applyAnnualHarvest(prePlan, prePlan, postPlan, one(c(1, NA, NA, NA)), elapsed = 0)
  expect_equal(terra::values(res$inventory)[, 2], c(0, 90, 100, 110))
})

test_that("no planned harvest leaves the inventory unchanged", {
  current <- inv(c(11, 11, 12, 12), c(81, 91, 101, 111))
  res <- applyAnnualHarvest(current, prePlan, postPlan, one(rep(NA, 4)), elapsed = 1)
  expect_equal(terra::values(res$inventory), terra::values(current))
  expect_equal(res$applied, 0)
})

test_that("harvestTiming defaults to 'period' (no change for existing runs)", {
  src <- paste(readLines(test_path("..", "..", "spades_ws3.R")), collapse = "\n")
  expect_match(src, 'defineParameter\\("harvestTiming", "character", "period"')
  expect_match(src, 'if \\(annual\\) sim <- applyAnnualHarvestSlice')
})
