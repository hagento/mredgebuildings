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

      if (grepl("A", tail(subSplit, 1))) {
        vars[["yStart"]] <- subSplit[[length(subSplit) - 2]]

        yEnd <- as.numeric(vars[["yStart"]]) + 4 # 5 year period
        vars[["yEnd"]] <- as.character(yEnd)

        vars[["subtype"]] <- sub("_A.nc", ".nc", subtype)
      }

      else if (grepl("B", tail(subSplit, 1))) {
        vars[["yEnd"]]   <- subSplit[[length(subSplit) - 1]]

        yStart <- as.numeric(vars[["yEnd"]]) - 4 # 5 year period
        vars[["yStart"]] <- as.character(yStart)

        vars[["yStartFile"]] <- subSplit[[length(subSplit) - 2]]

        vars[["subtype"]] <- sub("_B.nc", ".nc", subtype)
      }
    }

    else {stop("Invalid subtype given.")}

    return(vars)
  }


  # determine period range of subset
  getIdxRange <- function(yStart, yEnd, vars) {
    dStart <- as.Date(paste0(yStart, "-01-01"))
    dEnd   <- as.Date(paste0(yEnd,   "-12-31"))
    dRange <- seq.Date(from = dStart, to = dEnd, by = "day")

    if (is.null(vars[["yStartFile"]])) {
      idx <- seq(1:length(dRange))
    }
    else {
      dStartFile <- as.Date(paste0(vars[["yStartFile"]], "-01-01"))
      dRangeFile <- seq.Date(from = dStartFile, to = dEnd, by = "day")
      idx <- seq(length(dRangeFile) - length(dRange) + 1, length(dRangeFile))
    }

    return(idx)
  }


  # PROCESS DATA----------------------------------------------------------------

  vars <- splitSubtype(subtype)

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
    names(r) <- paste0("y", years[1]:years[2])

    # filter relevant years
    r <- terra::subset(r, as.numeric(substr(names(r), 2, 5)) > firstHistYear)

    # aggregate to common resolution of 0.5 deg
    if (any(raster::res(r) != 0.5)) {
      r <- terra::aggregate(r, fun = "sum",
                             fact = round(0.5 / terra::res(r), 3))
      terra::res(r) <- 0.5
      terra::ext(r) <- round(terra::ext(r))
    }

    x <- list(x = r, class = "SpatRaster", cache = FALSE)
  }


  else if (any(vars[["variable"]] %in% baitVars)) {

    if (!is.null(vars[["yStart"]])) {
      fpath <- file.path(vars[["variable"]], vars[["scenario"]], vars[["model"]], vars[["subtype"]])
      idx <- getIdxRange(vars[["yStart"]], vars[["yEnd"]], vars)
      r <- suppressWarnings(terra::rast(fpath, lyrs = idx))
    }
    else {
      fpath <- file.path(vars[["variable"]], vars[["scenario"]], vars[["model"]], subtype)
      r <- suppressWarnings(terra::rast(fpath))
    }

    x <- list(x = r, class = "SpatRaster", cache = FALSE)
  }

  else {stop("Subtype was incorrectly split or invalid subtype given.")}

  return(x)
}
