## spadesws3_params.py used to set `horizon = 10` unconditionally, after R had passed the
## spades_ws3 'horizon' parameter as py$horizon, so ws3 always planned 10 periods.

paramsFile <- test_path("..", "..", "python", "spadesws3_params.py")

horizonBlock <- function() {
  src <- readLines(paramsFile)
  i <- grep("^if 'horizon' not in globals", src)
  j <- grep("^horizon = int\\(horizon\\)", src)
  expect_length(i, 1); expect_length(j, 1)
  src[i:j]
}

test_that("spadesws3_params.py has no unconditional horizon assignment", {
  src <- readLines(paramsFile)
  expect_length(grep("^horizon\\s*=\\s*[0-9]", src), 0)
})

test_that("the horizon block keeps a value set from R and defaults to 10", {
  py <- Sys.which("python3")
  skip_if(!nzchar(py), "python3 not available")
  block <- horizonBlock()
  run <- function(pre) {
    f <- withr::local_tempfile(fileext = ".py")
    writeLines(c(pre, block, "print(horizon, type(horizon).__name__)"), f)
    system2(py, c("-I", f), stdout = TRUE)
  }
  expect_identical(run("horizon = 20.0"), "20 int")   # R numerics arrive as floats
  expect_identical(run("horizon = None"), "10 int")
  expect_identical(run(character(0)), "10 int")
})
