#' Read relevant ISIMIP data for mredgebuildings
#' 
#' @param subtype filename
#'
#' @author Hagen Tockhorn
#'
#' @importFrom stringr str_split
#'
#' NOTE:
#' folder structure in inputdata/sources/ISIMIPbuildings is expected to be:
#'    country masks : var/
#'    population :    var/scenario
#'    other :         var/scenario/model
#'    
#' NOTE: currently, this function only reads data from ISIMIP3b


readISIMIPbuildings <- function(subtype) {
  
  # PARAMETERS------------------------------------------------------------------
  
  baitVars <- c("tas", "sfcwind", "rsds", "huss")
  
  
  # FUNCTIONS-------------------------------------------------------------------
  
  splitSubtype <- function(subtype) {
    vars <- list()
    
    if ("countrymask" %in% subtype) {
      vars[["variable"]] = "countrymask"
    }
    
    else if (grepl("population", subtype)) {
      subSplit <- str_split(subtype, "_") %>% unlist()
      
      vars[["variable"]] <- subSplit[[1]]
      vars[["scenario"]] <- subSplit[[2]]
    }
    
    else if (any(baitVars %in% subtype)) {
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
  cwd <- getwd()
  
  if (vars[["variable"]] == "countrymask"){
    setwd("countrymask")
    varNames <- names(ncdf4::nc_open(subtype)[["var"]])
    countries <- list()
    for (var in varNames) {
      countries[[var]] <- raster::brick(files, varname = var)
    }
    r <- raster::brick(countries)
    names(r) <- gsub("m_", "", vars)
    setwd(cwd)
    x <- list(x = r, class = "RasterBrick")
  }
  
  
  if (vars[["variable"]] == "population") {
    setwd(file.path(vars[["variable"]], vars[["scenario"]]))
    r <- raster::brick(subtype)
    
    subtype <- gsub(".nc4", "", subtype)
    
    # rename years
    years <- strsplit(tail(strsplit(subtype, "_")[[1]], 1), "-")[[1]]
    names(r) <- paste0("y", years[1]:years[2])
    
    # filter relevant years
    r <- raster::dropLayer(
      r, as.numeric(substr(names(r), 2, 5)) < firstHistYear)
    
    # aggregate to common resolution of 0.5 deg
    if (any(raster::res(r) != 0.5)) {
      r <- raster::aggregate(r, fun = "sum",
                             fact = round(0.5 / raster::res(r), 3))
      raster::res(r)    <- 0.5
      raster::extent(r) <- round(raster::extent(r))
    }
    
    setwd(cwd)
    x <- list(x = raster::brick(r), class = "RasterBrick")
  }
  
  
  if (any(vars[["variable"]] %in% baitVars)) {
    setwd(file.path(vars[["variable"]], vars[["scenario"]], vars[["model"]]))
    r <- raster::brick(subtype, src = "") %>%
      round(digits=1)
    setwd(cwd)
    x <- list(x = r, class = "RasterBrick")
  }
  
  else {stop("Subtype was incorrectly split or invalid subtype given.")}
  
  return(x)
}
