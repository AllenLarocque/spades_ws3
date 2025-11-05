defineModule(sim, list(
  name = "spades_ws3",
  description = paste("This is the core module for spades_WS3 module family. It is a wrapper for WS3"),
  keywords = c("harvesting","dataInit","WS3"),
  authors = c(person("Gregory","Paradis", email = "0@01101.io", role = c("aut", "cre"))),
  childModules = character(0),
  version = list(SpaDES.core = "0.2.5.9000", spades_ws3 = "0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "spades_ws3.Rmd"),
  reqdPkgs = list('R.utils', 'reticulate'),
  parameters = rbind(
    defineParameter("mgmt.unit.theme", "numeric", NULL, NA, NA, "management unit theme index"),
    defineParameter("workers", "numeric", 1L, NA, NA, "number of worker threads to use for paralellized functions"),
    defineParameter("verbose", "numeric", 0, NA, NA, "console output verbosity level"),
    defineParameter("basenames", "character", NA, NA, NA, "MU baseneames to load"),
    defineParameter("enable.debugpy", "logical", FALSE, NA, NA, "enable debugpy"),
    defineParameter("horizon", "numeric", 1L, NA, NA, "ws3 simulation horizon (periods)"),
    defineParameter("base.year", "numeric", 2015L, NA, NA, "ws3 simulation base year"),
    defineParameter("scheduler.mode", "character", "optimize", NA, NA, "Switch between 'optimize' and 'areacontrol' harvest scheduler modes"),
    defineParameter("target.scalefactors", "numeric", NULL, NA, NA, "Target areas scale factors.  Only applicable if using 'areacontrol' scheduler mode."),
    defineParameter("mask.area.thresh", "numeric", 0., NA, NA, "Mask area threshold (for aggregation of bootstrapped masks).  Only applicable if using 'areacontrol' scheduler mode."),
    defineParameter("tif.path", 'character', 'tif', NA, NA, desc = 'name of directory with tifs in inputs'),
    defineParameter("yearOfFirstHarvest", 'numeric', start(sim), NA, NA, "year to schedule first harvest"),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between save events"),
    defineParameter(".useCache", "logical", FALSE, NA, NA, "Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "landscape", objectClass = "SpatRaster", desc = "stand age", sourceURL = NA)
  ),
  outputObjects = bind_rows(
    createsOutput(objectName = 'landscape', objectClass = 'SpatRaster', desc = 'SpatRaster landscape attributes')
  )
))

## event types

doEvent.spades_ws3 = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      sim <- Init(sim)
      sim <- scheduleEvent(sim, start(sim), "spades_ws3", "harvest")
      sim <- scheduleEvent(sim, start(sim), "spades_ws3", "grow")
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "spades_ws3", "plot")
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "spades_ws3", "save")
    },
    plot = {},
    save = {
      sim <- Save(sim)
      sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, "spades_ws3", "save")
    },
    harvest = {
      sim <- applyHarvest(sim)
      sim <- scheduleEvent(sim, time(sim) + 1, "spades_ws3", "harvest")
    },
    grow = {
      sim <- applyGrow(sim)
      sim <- scheduleEvent(sim, time(sim) + 1, "spades_ws3", "grow")
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

## event functions

Init <- function(sim) {

  if (is.null(P(sim)$basenames)) stop(paste("'basenames' parameter value not specified in", currentModule(sim)))
  cmp <- grep(pattern = paste0(currentModule(sim), "$"), x = list.files(modulePath(sim))) %>%
    list.files(path = modulePath(sim), full.names = TRUE)[.] # current module path

  py$dat_path<-inputPath(sim)
  py$basenames <- P(sim)$basenames
  py$enable_debugpy <- P(sim)$enable.debugpy
  py_run_file(file.path(cmp, "python", "spadesws3_params.py"))
  py$base_year <- P(sim)$base.year
  py$horizon <- P(sim)$horizon
  sim$fm <- py$bootstrap_forestmodel_kwargs()
  py$fm <- sim$fm
  return(invisible(sim))
}


Save <- function(sim) {
  sim <- saveFiles(sim)
  return(invisible(sim))
}


plotFun <- function(sim) {
  return(invisible(sim))
}


updateAges <- function(sim, offset = 1) {
  year <- as.integer(time(sim) - start(sim) + P(sim)$base.year)
  files1 <- sapply(P(sim)$basenames,
                   function(bn) file.path(inputPath(sim),
                                          P(sim)$tif.path,
                                          bn,
                                          paste("inventory_", toString(year), ".tif", sep="")))
  files2 <- sapply(P(sim)$basenames,
                   function(bn) file.path(inputPath(sim),
                                          P(sim)$tif.path,
                                          bn,
                                          paste("inventory_", toString(year+offset), ".tif", sep="")))


  rs.list <- lapply(files1, function(f) {
    r <- terra::rast(f)
    r <- terra::deepcopy(r)  # ensures a memory copy, not linked to disk
    r
  })
  names(rs.list) <- P(sim)$basenames  # Rename the list members their respective TSA names



  ###############################################################################
  # age.offset <- -1 # hack (why are age values in landscape raster stack off by 1?)
  ###############################################################################

  rs.list <- rapply(rs.list,  function(rs) {
    rs[[2]] <- terra::crop(sim$landscape$age, rs[[2]]) %>% terra::mask(., rs[[2]])
    rs[[2]][is.nan(rs[[2]])] <- NA
    return(rs)})
  mapply(terra::writeRaster, rs.list, files2, filetype='GTiff', overwrite=TRUE, datatype='INT4S')
  return(invisible(sim))
}


loadAges <- function(sim) {

  year <- as.integer(time(sim) - start(sim) + P(sim)$base.year)
  files <- sapply(P(sim)$basenames,
                  function(bn) file.path(inputPath(sim),
                                         P(sim)$tif.path,
                                         bn,
                                         paste("inventory_", toString(year), ".tif", sep="")))

  mergeAgeRasters <- function(files) {
    # Read band 2 from each raster
    rasters <- lapply(files, function(f) terra::rast(f, lyrs = 2))  # Age is in the second band

    if (length(rasters) > 1) {
      #r <- do.call(terra::mosaic, c(rasters, fun = "mean"))  # Merge by using `mosaic`, which is slower but handles overlapping cells
      r <- do.call(terra::merge, rasters) # merge by using `merge`, which is faster but may break with overlapping cells
    } else {
      r <- rasters[[1]]
    }

    # Replace NaN with NA
    # vals <- terra::values(r)
    # vals[is.nan(vals)] <- NA
    # terra::values(r) <- valsQ
    return(r)
  }
  r<-mergeAgeRasters(files)
  return(r)

}

applyHarvest <- function(sim) {
  year <- as.integer(time(sim) - start(sim) + P(sim)$base.year)
  py$base_year <- year
  sim$fm$base_year <- year
  updateAges(sim)
  py$simulate_harvest(fm = sim$fm,
                      basenames = P(sim)$basenames,
                      year = year,
                      mode = P(sim)$scheduler.mode,
                      target_scalefactors = P(sim)$target.scalefactors,
                      mask_area_thresh = P(sim)$mask.area.thresh,
                      verbose = P(sim)$verbose,
                      mgmt_unit_theme = P(sim)$mgmt.unit.theme,
                      workers=P(sim)$workers)

  sim$landscape[["age"]] <- loadAges(sim)

  return(invisible(sim))

}


applyGrow <- function(sim) {
  sim$landscape$age <- sim$landscape$age + 1
  updateAges(sim, offset=1)
  return(invisible(sim))
}


.inputObjects <- function(sim) {

  return(invisible(sim))
}
