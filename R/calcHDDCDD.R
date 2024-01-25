#' calculate HDD and CDD based on outdoor/indoor temperature difference
#'
#' @description heating and cooling degree days based on raw outside temperature
#'   or bais-adjusted internal temperature (BAIT), driver for space heating and
#'   cooling demand in buildings
#'
#' @param bait specify use of raw temperature or BAIT
#'
#' @return magpie object of heating and cooling degree days
#'
#' @author Robin Krekeler, Hagen Tockhorn
#' @examples
#'
#' \dontrun{
#' calcHDDCDD()
#' }
#'
#' @importFrom madrat getConfig
#' @importFrom raster brick cellStats subset stackApply getValues ncell
#'   xyFromCell res aggregate nlayers beginCluster endCluster
#' @importFrom ncdf4 nc_open
#' @importFrom tidyr %>%
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#' @importFrom pracma integral2
#' @importFrom stringr str_sub
#' @importFrom parallel mclapply


calcHDDCDD <- function(bait=FALSE) {
  
  # FILE MAPPING----------------------------------------------------------------
  
  # mappingFile <- "single-gcm_test.csv"
  
  mappingFile <- "cluster_test.csv"
  
  
  
  # COMPUTATION SETTINGS--------------------------------------------------------
  
  # number of cores to be considered in parallel computing
  ncores <- 16

  #raster::rasterOptions(tmpdir = "/p/tmp/hagento/rastertmp")
  

  # FUNCTIONS-------------------------------------------------------------------

  # fill dates for unnamed data
  #   ---> was supposed to go to convertISIMIP, doesn't work with country check
  fillDates <- function(r, filename, bait=FALSE) {
    if (grepl(".nc|.nc4", filename)) {
      filename <- gsub(".nc|.nc4", "", filename)
    }

    if (bait) {
      dStart <- as.Date(stringr::str_sub(filename, -17, -10),
                        format = "%Y%m%d")
      n <- raster::nlayers(r)
      dates <- seq.Date(dStart, by = "day", length.out = n)
    }
    else {
      yStart <- stringr::str_sub(filename, -9, -6)
      dStart <- as.Date(paste0(yStart, "-1-1"))
      n <- raster::nlayers(r)

      dates <- seq.Date(dStart, by = "day", length.out = n)
    }

    # fill dates
    names(r) <- dates
    return(r)
  }


  # prepare input for calcBAIT()
  calcBaitInput <- function(frsds=NULL, fsfc=NULL, fhuss=NULL, baitInput=NULL, mean=FALSE) {
    if(mean) {
      # optional: calculate daily means over years to fill missing data
      baitInputMean <- sapply(
        names(baitInput), function(var) {
          mean <- raster::stackApply(baitInput[[var]],
                                     unique(substr(names(baitInput[[var]]), 7, 11)),
                                     fun = "mean")
          names(mean) <- substr(names(mean), 7, 11)
          return(mean)})
    }
    else {
      input <- list(
        "rsds" = readSource("ISIMIPbuildings", subtype = frsds, convert = TRUE) %>%
          fillDates(frsds),
        "sfc"  = readSource("ISIMIPbuildings", subtype = fsfc, convert = TRUE) %>%
          fillDates(fsfc),
        "huss" = readSource("ISIMIPbuildings", subtype = fhuss, convert = TRUE) %>%
          fillDates(fhuss))
      return(input)
    }
  }


  #check if time period is congruent and adapt if necessary
  checkDates <- function(temp, baitInput) {
    dates_t <- substr(names(temp), 2, 11)

    baitInput <- sapply(names(baitInput), function(var) {
      # fill missing data with means from previous years
      # NOTE: "temp" and "baitInput" have the same global temporal lower
      #       boundary, since "temp" is the constraining dataset, only
      #       "baitInput" needs to be filled up.

      tmp <- baitInput[[var]]

      dates_b <- substr(names(tmp), 2, 11)

      datesFill <- setdiff(dates_t, dates_b)        # dates to fill up
      daysFill  <- unique(substr(datesFill, 6, 11))

      datesKeep <- intersect(dates_b, dates_t)      # dates to keep
      keep <- ifelse(length(datesKeep) > 0, TRUE, FALSE)

      if (keep) {
        tmp <- raster::subset(tmp, paste0("X", datesKeep))
        names(tmp) <- paste0("X", datesKeep)
        }

      if (length(daysFill) > 0) {
        baitInputMean <- calcBaitInput(mean = TRUE, baitInput = baitInput)

        # fill up missing dates with yearly-average value for specific day/cell
        baitInputFill <- raster::brick(lapply(
          daysFill, function(d){
            print(d)
            idx <- which(grepl(d, stringr::str_sub(datesFill, -5, -1)))
            r <- raster::brick(replicate(length(idx),
                                         baitInputMean[[var]][[paste0("X",d)]]))
            names(r) <- datesFill[idx]
            return(r)
            }
          )
        )

        # concatenate data
        if (keep) {tmp <- raster::brick(list(tmp, baitInputFill))}
        else      {tmp <- baitInputFill}

        # re-order dates
        tmp <- raster::brick(tmp[[order(names(tmp))]])
      }


      if (!identical(names(tmp), names(temp))) {
        return(print("Warning: Dates of Temperature and BAIT Input Data are not aligned."))
      }
      return(tmp)
      },
      USE.NAMES = TRUE
    )
    return(baitInput)
  }


  #--- CALCULATE BAIT (building-adjusted internal temperature)

  # NOTE: The parameters are taken from Staffell et al. 2023

  # counterfactuals for solar, wind, humidity, temperature
  cfac <- function(t, type, params=NULL) {
    if (is.null(params)) {
      params <- switch(type,
                       s = c(100, 7),
                       w = c(4.5, 0.025),
                       h = c(1.1, 0.06),
                       t = c(16))}

    if      (type == "s") {return(params[[1]] + params[[2]]*t)}
    else if (type == "w") {return(params[[1]] - params[[2]]*t)}
    else if (type == "h") {return(exp(params[[1]] + params[[2]]*t))}
    else if (type == "t") {return(params[[1]])}

    else {print("No valid parameter type specified.")}
  }


  calcBAIT <- function(temp, baitInput, weight=NULL, params=NULL) {
    # weight=(x,y,z,sigma,bLower,bUpper,bMax)
    if (is.null(weight)) {
      print("Please give appropriate weights for the calculation of BAIT.")
      weight <- c(1,1,1,1,1,1,1)}

    dates <- names(temp)

    solar <- baitInput$rsds %>% brick()
    wind  <- baitInput$sfc  %>% brick()
    hum   <- baitInput$huss %>% brick()

    s <- solar -  cfac(temp, type="s", params)
    w <- wind  -  cfac(temp, type="w", params)
    h <- hum   -  cfac(temp, type="h", params)
    t <- temp  -  cfac(temp, type="t", params)

    # calc raw bait
    bait <- temp + weight[[1]]*s - weight[[2]]*w + weight[[3]]*h*t

    # smooth bait over preceding two days with smoothing parameter sigma
    bait <- raster::brick(
      sapply(seq(raster::nlayers(bait)), function(i) {
        if (i %in% c(1, 2)) {return(bait[[i]])}
        else {
          bait[[i]] <- (bait[[i]] + weight[[4]]*bait[[i-1]] + weight[[4]]**2 * bait[[i-2]]) / (1 + weight[[4]] + weight[[4]]**2)
          return(bait[[i]])
          }
        }
      )
    )

    # weighted blend of BAIT and raw temperature
    bBar <- (temp - 0.5*(weight[[6]] + weight[[5]])) * 10 / (weight[[6]] - weight[[5]])
    b <- weight[[7]] / (1 + exp(-bBar))

    bait <- bait * (1 - b) + (temp * b)

    names(bait) <- dates
    return(bait)
  }


  #--- CALCULATE HDD/CDD

  # Calculate heating factor matrix for normally distributed ambient
  # and limit temperatures

  # t1 : ambient temperature variable
  # t2 : limit temperature variable

  heatingFactor <- function(t2, t1, tamb, tamb_std, tlim, tlim_std) {
    h <- dnorm(t2, mean=tlim, sd=tlim_std) * dnorm(t1, mean=tamb, sd=tamb_std) * (t2 - t1)
    return(h)
  }

  coolingFactor <- function(t2, t1, tamb, tamb_std, tlim, tlim_std) {
    h <- dnorm(t2, mean=tlim, sd=tlim_std) * dnorm(t1, mean=tamb, sd=tamb_std) * (t1 - t2)
    return(h)
  }

  # check if ambient/limit temperature interval is reasonable
  # e.g. t_lim = 17C and t_amb = -50C wouldn't give reasonable CDD
  checkTDif <- function(tamb, tlim, typeDD, tamb_std, tlim_std) {
    check <- TRUE
    stdDif <- tamb_std + tlim_std
    if (typeDD == "HDD") {
      if (tamb - tlim > stdDif) {
        check <- FALSE
      }
    }
    else if (typeDD == "CDD") {
      if (tlim - tamb > 2*stdDif) {
        check <- FALSE
      }
    }
    return(check)
  }

  # calculate HDD/CDD values per day for given ambient/limit temp combination
  calcHDDCDDFactors <- function(tlow, tup, tlim, tamb_std=5, tlim_std=5) {
    t <- seq(tlow, tup, .1)

    hddcddFactors <- do.call(
      "rbind", lapply(
        c("HDD", "CDD"), function(typeDD) {
          do.call(
            "rbind", lapply(
              t, function(tamb) {
                do.call(
                  "rbind", lapply(
                    tlim[[typeDD]], function(.tlim) {
                      if (!checkTDif(tamb, .tlim, typeDD, tamb_std, tlim_std)) {
                        tmp <- data.frame("T_amb"        = tamb,
                                          "T_amb_K"      = round(tamb + 273.15, 1),
                                          "T_lim"        = .tlim,
                                          "factor"       = 0,
                                          "factor_err"   = 0,
                                          "typeDD"       = typeDD)
                      }
                      else {

                        # integration boundaries
                        x1 <- .tlim - 4*tlim_std
                        x2 <- .tlim + 4*tlim_std
                        y1 <- min(.tlim - 3*tlim_std, tamb - 3*tlim_std)
                        y2 <- max(.tlim + 3*tlim_std, tamb + 3*tlim_std)

                        if (typeDD == "HDD") {
                          f <- pracma::integral2(heatingFactor,
                                                 xmin = x1,
                                                 xmax = x2,
                                                 ymin = y1,
                                                 ymax = function(x){x},
                                                 tamb = tamb,
                                                 tamb_std = tamb_std,
                                                 tlim = .tlim,
                                                 tlim_std = tlim_std,
                                                 reltol = 1e-1
                          )
                        }
                        else {
                          f <- pracma::integral2(coolingFactor,
                                                 xmin = x1,
                                                 xmax = x2,
                                                 ymin = function(x){x},
                                                 ymax = y2,
                                                 tamb = tamb,
                                                 tamb_std = tamb_std,
                                                 tlim = .tlim,
                                                 tlim_std = tlim_std,
                                                 reltol = 1e-1)
                        }
                        tmp <- data.frame("T_amb"        = tamb,
                                          "T_amb_K"      = round(tamb + 273.15, 1),
                                          "T_lim"        = .tlim,
                                          "factor"       = f$Q,
                                          "factor_err"   = f$error,
                                          "typeDD"       = typeDD)
                      }
                    }
                  )
                )
              }
            )
          )
        }
      )
    )
    return(hddcddFactors)
  }



  # fill HDD/CDD from factors for given ambient/limit temperature combination
  calcCellHDDCDD <- function(.temp, .typeDD, .tlim, factors) {
    # extract years
    years <- names(.temp) %>%
      substr(2, 5) %>%
      as.numeric()


    # add tolerance of 0.04K to avoid machine precision errors
    factors <- factors %>%
      filter(.data[["typeDD"]] == .typeDD, .data[["T_lim"]] == .tlim) %>%
      dplyr::reframe(from = .data[["T_amb_K"]] - 0.04,
                     to = .data[["T_amb_K"]] + 0.04,
                     becomes = .data[["factor"]]) %>%
      data.matrix()


    # swap ambient temperature values with corresponding DD values
    hddcdd <- raster::reclassify(.temp, factors)

    # aggregate to yearly HDD/CDD [K.d/a]
    hddcdd <- raster::stackApply(hddcdd, years, fun = sum)
    names(hddcdd) <- gsub("index_", "y", names(hddcdd))

    return(hddcdd)
  }


  # aggregate cellular HDD/CDD to country-wide average (population-weighted)
  aggCells <- function(r, weight, mask) {
    years_r <- names(r)
    years_w <- names(weight)

    if (!all(years_r %in% years_w)) {
      stop("Time periods of raster file and aggregation weights do not match.")
    }

    # loop: years in raster file r
    hddcdd_agg <- do.call(
      "rbind", lapply(
        years_r, function(y) {
          tmp <- raster::subset(r, y) * raster::subset(weight, y) * mask
          tmp_tot <- raster::subset(weight, y) * mask
          tmp <- raster::cellStats(tmp, "sum") /
            raster::cellStats(tmp_tot, "sum")
          tmp <- data.frame("region" = names(mask),
                            "period" = y,
                            "value"  = round(tmp, 1))
          rownames(tmp) <- c()
          return(tmp)
        }
      )
    )
  }


  # calculate all desired output from given tas file
  calcStackHDDCDD <- function(file, tlim, countries, pop, factors, bait,
                              frsds = NULL,
                              fsfc  = NULL,
                              fhuss = NULL,
                              wBAIT = NULL) {
   
    #raster::removeTmpFiles()
 
    # read cellular temperature
    temp <- readSource("ISIMIPbuildings", subtype = file, convert = TRUE) %>%
      fillDates(file)

    dates <- names(temp)

    # optional: transform raw temperature into BAIT
    if (bait) {
      baitInput <- calcBaitInput(frsds, fsfc, fhuss)

      # note: easier to do in [C] then convert back
      temp <- temp - 273.15   # [C]
      print("check dated bait")
      baitInput <- checkDates(temp, baitInput)
      print("calc bait")
      temp <- calcBAIT(temp, baitInput, weight = wBAIT)
      temp <- round(temp + 273.15, digits=1)    # [K]
    }

    names(temp) <- dates

    # loop: typeDD
    hddcdd <- do.call(
      "rbind", lapply(
        c("HDD", "CDD"), function(typeDD) {
          # loop: threshold temperatures
          do.call(
            "rbind", lapply(
              tlim[[typeDD]], function(t) {
                hddcdd_agg <- calcCellHDDCDD(temp, typeDD, t, factors) %>%
                  aggCells(pop, countries) %>%
                  mutate("variable" = typeDD,
                         "tlim"     = t)    # [C]
              }
            )
          )
        }
      )
    )
  }


  plotRaster <- function(r) {
    library("ggplot2")
    r_map <- stack(as.data.frame(raster::getValues(r)))
    coords <- raster::xyFromCell(r, seq_len(raster::ncell(r)))
    names(r_map) <- c('value', 'variable')
    r_map <- cbind(coords, r_map)
    ggplot(r_map) +
      geom_tile(aes(x, y, fill = value)) +
      facet_wrap(~ variable) +
      scale_fill_gradientn(colours = rev(terrain.colors(225))) +
      coord_equal()
  }

  roundRas <- function(a) {
    a <- round(a, digits = 1) %>%
      brick()
  }



  # PARAMETERS------------------------------------------------------------------

  # threshold temperature for heating and cooling [C]
  # NOTE: Staffel gives global average of T_heat = 14, T_cool = 20
  # t_lim <- list("HDD" = seq(17, 25), "CDD" = seq(17, 25))
  t_lim <- list("HDD" = seq(17, 17), "CDD" = seq(23, 23))

  # standard deviations for temperature distributions
  tlim_std <- 5   # threshold
  tamb_std <- 5   # ambient

  # historical years    -> NOTE: not sure if further needed
  firstHistYear <- 1950
  lastHistYear  <- 2009

  # range of pre-calculated HDD/CDD-values, e.g. [223, 232] K, converted to [C]
  tlow <- 223 - 273.15
  tup <- 323 - 273.15



  #--- BAIT parameters

  # The weights (x,y,z) for calcBAIT and the smoothing coefficient are assumed to
  # be region-independent and equal to the mean of the values given in Staffell
  # et al. 2023
  x <- 0.012
  y <- -0.20
  z <- 0.05

  # smoothing coefficient
  sig <- 0.50

  # The blending parameters for the blending of BAIT and raw temperature are like-wise
  # taken from the paper.
  bLower <- 15
  bUpper <- 23
  bMax   <- 0.5

  # concatenate to vector
  wBAIT <- c(x, y, z, sig, bLower, bUpper, bMax)



  # READ-IN DATA----------------------------------------------------------------
  
  # list of files that are processed
  files <- toolGetMapping(mappingFile, type = "sectoral", where = "mappingfolder") %>%
    filter(variable != "")
  
  # cells -> country
  fCM <- file.path(files[files$variable == "CountryMask", "file"])
  countries <- readSource("ISIMIPbuildings", subtype = fCM, convert = FALSE)



  # PROCESS DATA----------------------------------------------------------------

  # calculate HDD/CDD-factors
  hddcddFactor <- calcHDDCDDFactors(tlow=-50.15, tup=49.85, t_lim, tamb_std, tlim_std)

  # loop: GCM results for ambient temperature (SSP scenarios)
  hddcdd_multi <- do.call(
    "rbind",
    mclapply( # ssp iteration
      ssps <-  files %>%
        filter(.data[["variable"]] == "tas") %>%
        select("ssp") %>%
        unique() %>%
        as.list(),
      function(s) {
        fpop <- files %>% filter(ssp == s, variable == "pop")
        pop <- readSource("ISIMIPbuildings", subtype = fpop$file,
                          convert = FALSE)
        do.call(
          "rbind",
          lapply( # rcp iteration
            rcps <- files %>%
              filter(.data[["variable"]] == "tas",
                     .data[["ssp"]] == s) %>%
              select("rcp") %>%
              unique() %>%
              as.list(),
            function(r) {
              do.call(
                "rbind",
                lapply( # model iteration
                  model <- files %>%
                    filter(.data[["variable"]] == "tas") %>%
                    select("gcm") %>%
                    unique() %>%
                    as.list(),
                  function(m) {
                    f <- filter(files, .data[["ssp"]] == s, .data[["rcp"]] == r)
                    do.call( # file iteration
                      "rbind",
                      lapply(
                        seq(nrow(filter(f, f$variable == "tas"))),
                        function(n) {
                          ftas <- file.path(f[f$variable == "tas" & f$gcm == m,][[n, "file"]])
                          
                          print(paste("Processing temperature file:", ftas))

                          if(bait) {
                            frsds <- file.path(f[f$variable == "rsds",][[n, "file"]])
                            
                            fsfc  <- file.path(f[f$variable == "sfc",][[n, "file"]])
                            
                            fhuss <- file.path(f[f$variable == "huss",][[n, "file"]])

                            hddcddCell <- calcStackHDDCDD(ftas,
                                                          t_lim,
                                                          countries,
                                                          pop,
                                                          hddcddFactor,
                                                          bait,
                                                          frsds = frsds,
                                                          fsfc  = fsfc,
                                                          fhuss = fhuss,
                                                          wBAIT = wBAIT)}

                          else {
                            hddcddCell <- calcStackHDDCDD(ftas,
                                                          t_lim,
                                                          countries,
                                                          pop,
                                                          hddcddFactor,
                                                          bait)}

                          hddcddCell <- hddcddCell %>%
                            mutate("model" = m,
                                   "ssp" = s,
                                   "rcp" = r)
                        }
                      )
                    )
                  }
                )
              )
            }
          )
        )
      }, 
      mc.cores = ncores
    )
  )

  rownames(hddcdd) <- c()

  # average over all GCM's
  # data <- hddcdd %>%
  #   group_by(across(-all_of(c("model","value")))) %>%
  #   summarise(value = mean(.data[["value"]]))



  # OUTPUT----------------------------------------------------------------------

  # might be necessary for magpie output
  data <- hddcdd %>%
    mutate(period = gsub("y", "", .data[["period"]])) %>%
    unite("variable", .data[["variable"]], .data[["tlim"]], sep = "_") %>%
    unite("scenario", .data[["ssp"]], .data[["rcp"]], sep = "_")

  # save data as csv
  # optFolder <- getConfig("outputfolder")
  # write.csv(data, file = file.path(optFolder, "hddcdd.csv"), row.names = FALSE)


  return(list(x = data))
}



