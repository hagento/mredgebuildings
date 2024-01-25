#' Read relevant ISIMIP data for mredgebuildings
#' 
#' @param subtype filename
#'
#' @author Hagen Tockhorn
#'
#' @importFrom stringr str_split
#' @importFrom raster brick dropLayer res extent aggregate
#'
# NOTE:
# folder structure in inputdata/sources/ISIMIPbuildings is expected to be:
#    country masks : var/
#    population :    var/scenario
#    other :         var/scenario/model
#    
# NOTE: currently, this function only reads data from ISIMIP3b


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
    }
    
    else {stop("Invalid subtype given.")}
    
    return(vars)
  }
  
  
  # PROCESS DATA----------------------------------------------------------------
  
  vars <- splitSubtype(subtype)
  
  if (vars[["variable"]] == "countrymask"){
    fpath <- file.path("countrymasks", subtype)
    varNames <- names(ncdf4::nc_open(fpath)[["var"]])
    countries <- list()
    for (var in varNames) {
      countries[[var]] <- suppressWarnings(brick(fpath, varname = var))
    }
    r <- brick(countries)
    names(r) <- gsub("m_", "", varNames)
    
    x <- list(x = r, class = "RasterBrick")
  }
  
  
  else if (vars[["variable"]] == "population") {
    fpath <- file.path(vars[["variable"]], vars[["scenario"]], subtype)
    r <- suppressWarnings(brick(fpath))
    
    subtype <- gsub(".nc", "", subtype)
    
    # rename years
    years <- tail(strsplit(subtype, "_")[[1]], 2)
    names(r) <- paste0("y", years[1]:years[2])
    
    # filter relevant years
    r <- dropLayer(
      r, as.numeric(substr(names(r), 2, 5)) < firstHistYear)
    
    # aggregate to common resolution of 0.5 deg
    if (any(raster::res(r) != 0.5)) {
      r <- raster::aggregate(r, fun = "sum",
                             fact = round(0.5 / raster::res(r), 3))
      raster::res(r)    <- 0.5
      raster::extent(r) <- round(raster::extent(r))
    }
    
    x <- list(x = brick(r), class = "RasterBrick")
  }
  
  
  else if (any(vars[["variable"]] %in% baitVars)) {
    fpath <- file.path(vars[["variable"]], vars[["scenario"]], vars[["model"]], subtype)
    r <- suppressWarnings(brick(fpath, src = "")) %>%
      round(digits=1)
    
    x <- list(x = r, class = "RasterBrick")
  }
  
  else {stop("Subtype was incorrectly split or invalid subtype given.")}
  
  return(x)
}
