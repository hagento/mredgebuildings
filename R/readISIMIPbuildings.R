#' Read relevant ISIMIP data for mredgebuildings
#'
#'
#' @param subtype filename
#'
#' @author Hagen Tockhorn
#'
#' @importFrom stringr str_split
#' @importFrom raster brick dropLayer res extent aggregate
#'
#' @note
#' Argument subtype consist of full filenames of the respective raster files. In
#' order to improve runtime and general problems that come with processing large
#' raster files, filenames marked with "_A" or "_B" before the filetype declaration,
#' will be split in the respective first and second half of the dataset.
#'
#' @note
#' folder structure in inputdata/sources/ISIMIPbuildings is expected to be:
#'    country masks : var/
#'    population :    var/scenario
#'    other :         var/scenario/model
#'
#' @note currently, this function only reads data from ISIMIP3b


readISIMIPbuildings <- function(subtype) {

  # PARAMETERS------------------------------------------------------------------

  baitVars <- c("tas", "sfcwind", "rsds", "huss")

  firstHistYear <- 1960


  # FUNCTIONS-------------------------------------------------------------------

  splitSubtype <- function(subtype) {
    vars <- list()

    if (grepl("countrymask", subtype)) {
      vars[["variable"]] = "countrymask"
    }

    else if (grepl("population", subtype)) {
      subSplit <- str_split(subtype, "_") %>% unlist()

      vars[["variable"]] <- subSplit[[1]]
      vars[["scenario"]] <- subSplit[[2]]
    }

    else if (any(sapply(baitVars, grepl, x = subtype))) {
      subSplit <- str_split(subtype, "_") %>% unlist()

      vars[["variable"]] <- subSplit[[5]]
      vars[["scenario"]] <- subSplit[[4]]
      vars[["model"]]    <- subSplit[[1]]

      # raster data will be split into individual years
      if (length(subSplit) > 9) {
        # split index defines the year
        vars[["idx"]] <- gsub(".nc", "", tail(subSplit, 1)) %>%
          as.numeric()

        # temporal range of data
        vars[["yStart"]] <- subSplit[[8]]
        vars[["yEnd"]]   <- subSplit[[9]]

        # year of interest
        vars[["year"]] <- seq(vars[["yStart"]] %>%
                                as.numeric(),
                              vars[["yEnd"]] %>%
                                as.numeric())[[vars[["idx"]]]] %>%
          as.character()

        vars[["subtype"]] <- sub("_(1[0-9]|\\d)\\.nc$", ".nc", subtype)
      }
    }

    else {stop("Invalid subtype given.")}

    return(vars)
  }


  # determine period range of subset
  getRanges <- function(vars) {
    # total temporal range
    dRange <- seq.Date(from = as.Date(paste0(vars[["yStart"]], "-01-01")),
                       to = as.Date(paste0(vars[["yEnd"]],   "-12-31")),
                       by = "day")

    # temporal range of interest
    yRange <- seq.Date(from = as.Date(paste0(vars[["year"]], "-01-01")),
                       to   = as.Date(paste0(vars[["year"]], "-12-31")),
                       by   = "day")

    # indices of range of interest
    idxRange <- match(yRange, dRange)

    return(list("yRange"   = yRange,
                "idxRange" = idxRange))
  }


  # PROCESS DATA----------------------------------------------------------------
  vars <- splitSubtype(subtype)

  # region mask
  if (vars[["variable"]] == "countrymask"){
    fpath <- file.path("countrymasks", subtype)
    varNames <- names(ncdf4::nc_open(fpath)[["var"]])
    countries <- list()
    for (var in varNames) {
      countries[[var]] <- suppressWarnings(terra::rast(fpath, subds = var))
    }
    r <- terra::rast(countries)
    names(r) <- gsub("m_", "", varNames)

    x <- list(x = r, class = "SpatRaster")
  }


  # population
  else if (vars[["variable"]] == "population") {
    fpath <- file.path(vars[["variable"]], vars[["scenario"]], subtype)

    if (vars[["scenario"]] == "picontrol") {
      r <- suppressWarnings(terra::rast(fpath))
    }
    else {
      r <- suppressWarnings(terra::rast(fpath, subds = "total-population"))
    }

    subtype <- gsub(".nc", "", subtype)

    # rename years
    years <- tail(strsplit(subtype, "_")[[1]], 2)
    names(r) <- years[1]:years[2]

    # filter relevant years
    r <- terra::subset(r, as.numeric(names(r)) >= firstHistYear)

    # aggregate to common resolution of 0.5 deg
    if (any(raster::res(r) != 0.5)) {
      r <- terra::aggregate(r, fun = "sum",
                             fact = round(0.5 / terra::res(r), 3))
      terra::res(r) <- 0.5
      terra::ext(r) <- round(terra::ext(r))
    }

    x <- list(x = r, class = "SpatRaster", cache = FALSE)
  }


  # climate data
  else if (any(vars[["variable"]] %in% baitVars)) {
    # slice single years
    if (!is.null(vars[["yStart"]])) {
      fpath  <- file.path(vars[["variable"]], vars[["scenario"]], vars[["model"]], vars[["subtype"]])
      ranges <- getRanges(vars)

      r        <- suppressWarnings(terra::rast(fpath, lyrs = ranges[["idxRange"]]))
      names(r) <- ranges[["yRange"]]
    }

    # full data set
    else {
      fpath <- file.path(vars[["variable"]], vars[["scenario"]], vars[["model"]], subtype)
      r <- suppressWarnings(terra::rast(fpath))
    }

    x <- list(x = r, class = "SpatRaster", cache = FALSE)
  }

  else {stop("Subtype was incorrectly split or invalid subtype given.")}

  return(x)
}
