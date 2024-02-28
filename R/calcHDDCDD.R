#' calculate HDD and CDD based on outdoor/indoor temperature difference
#'
#' @description heating and cooling degree days based on raw outside temperature
#'   or bais-adjusted internal temperature (BAIT), driver for space heating and
#'   cooling demand in buildings
#'
#' @param mappingFile data.frame containing input data file names and directories
#' @param bait specify use of raw temperature or BAIT
#' @param multiscen specify if mappingFile covers more than one scenario
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
#' @importFrom raster cellStats
#' @importFrom ncdf4 nc_open
#' @importFrom tidyr %>%
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#' @importFrom pracma integral2
#' @importFrom stringr str_sub
#' @importFrom terra app nlyr tapp subset rast classify time


calcHDDCDD <- function(mappingFile, bait=FALSE, multiscen = FALSE) {


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
      n <- nlyr(r)
      dates <- seq.Date(dStart, by = "day", length.out = n)
    }
    else {
      yStart <- stringr::str_sub(filename, -9, -6)
      dStart <- as.Date(paste0(yStart, "-1-1"))
      n <- nlyr(r)

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
          mean <- tapp(baitInput[[var]],
                              unique(substr(names(baitInput[[var]]), 6, 11)),
                              fun = "mean")
          names(mean) <- gsub("\\.", "-", substr(names(mean), 2, 6))
          return(mean)})
    }
    else {
      input <- list(
        "rsds" = readSource("ISIMIPbuildings", subtype = frsds, convert = TRUE),
        "sfc"  = readSource("ISIMIPbuildings", subtype = fsfc,  convert = TRUE),
        "huss" = readSource("ISIMIPbuildings", subtype = fhuss, convert = TRUE))
      return(input)
    }
  }


  #check if time period is congruent and adapt if necessary
  checkDates <- function(baitInput, temp) {
    dates_t <- names(temp)

    baitInput <- sapply(names(baitInput), function(var) {
      # fill missing data with means from previous years
      # NOTE: "temp" and "baitInput" have the same global temporal lower
      #       boundary, since "temp" is the constraining dataset, only
      #       "baitInput" needs to be filled up.

      tmp <- baitInput[[var]]

      dates_b <- names(tmp)

      datesFill <- setdiff(dates_t, dates_b)        # dates to fill up
      daysFill  <- unique(substr(datesFill, 6, 11))

      datesKeep <- intersect(dates_b, dates_t)      # dates to keep
      keep <- ifelse(length(datesKeep) > 0, TRUE, FALSE)

      if (keep) {
        tmp <- subset(tmp, datesKeep)
        names(tmp) <- datesKeep
        }

      if (length(daysFill) > 0) {
        baitInputMean <- calcBaitInput(mean = TRUE, baitInput = baitInput)

        # fill up missing dates with yearly-average value for specific day/cell
        baitInputFill <- rast(lapply(
          daysFill, function(d){
            print(d)
            idx <- which(grepl(d, stringr::str_sub(datesFill, -5, -1)))
            r <- rast(replicate(length(idx),
                                         baitInputMean[[var]][[d]]))
            names(r) <- datesFill[idx]
            return(r)
            }
          )
        )

        # concatenate data
        if (keep) {tmp <- rast(list(tmp, baitInputFill))}
        else      {tmp <- baitInputFill}

        # re-order dates
        tmp <- rast(tmp[[order(names(tmp))]])
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
                       w = c(4.5, -0.025),
                       h = c(1.1, 0.06),
                       t = c(16))}

    print(names(params))
    if      (type == "s") {return(params[[1]] + params[[2]]*t)}
    else if (type == "w") {return(params[[1]] + params[[2]]*t)}
    else if (type == "h") {return(exp(params[[1]] + params[[2]]*t))}
    else if (type == "t") {return(params[[1]])}

    else {print("No valid parameter type specified.")}
  }


  # smooth data
  smooth <- function(r, weight) {
    # smooth bait over preceding two days with smoothing parameter sigma
    print("smooth")

    # one day indented
    r1D <- r[[c(nlyr(r), 1:(nlyr(r) - 1))]]
    r1D[[1]] <- 0

    # two days indented
    r2D <-  r[[c(nlyr(r)-1, nlyr(r), 1:(nlyr(r) - 2))]]
    r2D[[1:2]] <- 0

    # smooth
    rSmooth <- (r + weight[[4]]*r1D + weight[[4]]**2 * r2D) / (1 + weight[[4]] + weight[[4]]**2)

    return(rSmooth)
  }

  blend <- function(baitDF, temp, weight) {
    # weighted blend of BAIT and raw temperature
    print("blend")
    bBar <- (temp - 0.5*(weight[[6]] + weight[[5]])) * 10 / (weight[[6]] - weight[[5]])
    b <- weight[[7]] / (1 + exp(-bBar))

    baitDF <- baitDF * (1 - b) + (temp * b)
    return(baitDF)
  }


  calcBAIT <- function(baitInput, temp, weight=NULL, params=NULL) {
    # weight=(x,y,z,sigma,bLower,bUpper,bMax)
    if (is.null(weight)) {
      print("Please give appropriate weights for the calculation of BAIT.")
      weight <- c(1,1,1,1,1,1,1)}

    dates <- names(temp)

    solar <- baitInput$rsds
    wind  <- baitInput$sfc
    hum   <- baitInput$huss

    print("calc s")
    s <- solar -  cfac(temp, type="s", params = c(params[["a_rsds"]], params[["b_rsds"]]))
    print("calc w")
    w <- wind  -  cfac(temp, type="w", params = c(params[["a_sfcwind"]], params[["b_sfcwind"]]))
    print("print h")
    h <- hum   -  cfac(temp, type="h", params = c(params[["a_huss"]], params[["b_huss"]]))
    print("calc t")
    t <- temp  -  cfac(temp, type="t", params = NULL)

    # calc raw bait
    print("calc bait")
    bait <- temp + weight[[1]]*s + weight[[2]]*w + weight[[3]]*h*t

    bait <- smooth(bait, weight)
    bait <- blend(bait, temp, weight)

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
  calcCellHDDCDD <- function(temp, .typeDD, .tlim, factors) {
    # extract years
    dates <- names(temp)

    # add tolerance of 0.04K to avoid machine precision errors
    factors <- factors %>%
      filter(.data[["typeDD"]] == .typeDD, .data[["T_lim"]] == .tlim) %>%
      dplyr::reframe(from = .data[["T_amb_K"]] - 0.049,
                     to = .data[["T_amb_K"]] + 0.049,
                     becomes = .data[["factor"]]) %>%
      data.matrix()

    # swap ambient temperature values with corresponding DD values
    hddcdd <- classify(temp, factors)

    terra::time(hddcdd) <- as.Date(dates)

    # aggregate to yearly HDD/CDD [K.d/a]
    hddcdd <- tapp(hddcdd, "years", fun = sum, na.rm = TRUE)

    names(hddcdd) <- gsub("_", "", names(hddcdd))
    return(hddcdd)
  }


  # aggregate cellular HDD/CDD to country-wide average (population-weighted)
  aggCells <- function(r, weight, mask) {
    years_r <- names(r)
    years_w <- names(weight)

    print(paste("years_r:", years_r))
    print(paste("years_w:", years_w))

    if (!all(years_r %in% years_w)) {
      stop("Time periods of raster file and aggregation weights do not match.")
    }

    # loop: years in raster file r
    hddcddAgg <- do.call(
      "rbind", lapply(
        years_r, function(y) {
          browser()
          tmp <- subset(r, y) * subset(weight, y) * mask
          tmpTot <- subset(weight, y) * mask

          tmpSum <- terra::global(tmp, "sum", na.rm = TRUE)$sum
          tmpTotSum <- terra::global(tmpTot, "sum", na.rm = TRUE)$sum
          tmp <- tmpSum / tmpTotSum
          tmp <- data.frame("region" = names(mask),
                            "period" = y,
                            "value"  = round(tmp, 3))
          rownames(tmp) <- c()
          return(tmp)
        }
      )
    )
    return(hddcddAgg)
  }


  # calculate all desired output from given tas file
  calcStackHDDCDD <- function(file, tlim, countries, pop, factors, bait,
                              frsds = NULL,
                              fsfc  = NULL,
                              fhuss = NULL,
                              wBAIT = NULL,
                              params = NULL) {

    # read cellular temperature
    temp <- readSource("ISIMIPbuildings", subtype = file, convert = TRUE)

    dates <- names(temp)

    # optional: transform raw temperature into BAIT
    if (bait) {
      # note: easier to do in [C]
      temp <- temp - 273.15   # [C]

      # read and prepare bait input data
      baitInput <- calcBaitInput(frsds, fsfc, fhuss) %>%
        checkDates(temp)

      # calculate bait
      temp <- calcBAIT(baitInput, temp, weight = wBAIT, params = params)

      # convert back to [K]
      temp <- temp + 273.15   # [K]
    }

    temp <- terra::round(temp, digits = 1)
    names(temp) <- dates

    print("Calculating HDD/CDDs per cell.")

    fSplit <- str_split(file, "_") %>% unlist()

    # loop: typeDD
    hddcdd <- do.call(
      "rbind", lapply(
        c("HDD", "CDD"), function(typeDD) {
          # loop: threshold temperatures
          do.call(
            "rbind", lapply(
              tlim[[typeDD]], function(t) {
                hddcdd_agg <- calcCellHDDCDD(temp, typeDD, t, factors)

                # save intermediate results
                decades <- c(as.numeric(gsub("y", "", names(hddcdd_agg))) %% 10 == 0)
                if (any(decades)) {
                  hddcdd_save <- hddcdd_agg[[decades]]

                  y <- names(hddcdd_agg)[as.numeric(gsub("y", "", names(hddcdd_agg))) %% 10 == 0]

                  y <- gsub("y", "", y) %>%
                    paste(collapse = "-")

                  rname <- paste0(fSplit[[1]], "_", y, "_", fSplit[[4]], "_", typeDD, "_", t)

                  if (bait) {rname <- paste0(rname, "_bait.nc")}
                  else {rname <- paste0(rname, ".nc")}


                  terra::writeRaster(hddcdd_save,
                                     paste0("/p/tmp/hagento/output/rasterdata/", rname),
                                     overwrite = TRUE)
                }

                hddcdd_agg <- hddcdd_agg %>%
                  aggCells(pop, countries) %>%
                  mutate("variable" = typeDD,
                         "tlim"     = t)    # [C]

                return(hddcdd_agg)
              }
            )
          )
        }
      )
    )
  }


  # initialize full calculation
  makeCalculations <- function(f, m, n, t_lim, countries, pop, hddcddFactor,
                               bait, wBAIT = NULL, params = NULL, suffix = NULL) {

    if (!is.null(suffix)) {
      ftas <- gsub(".nc",
                   paste0("_", suffix, ".nc"),
                   f[f$variable == "tas"       & f$gcm == m,][[n, "file"]])

      if (bait){
        frsds <- gsub(".nc",
                      paste0("_", suffix, ".nc"),
                      f[f$variable == "rsds"    & f$gcm == m,][[n, "file"]])

        fsfc  <- gsub(".nc",
                      paste0("_", suffix, ".nc"),
                      f[f$variable == "sfcwind" & f$gcm == m,][[n, "file"]])

        fhuss <- gsub(".nc",
                      paste0("_", suffix, ".nc"),
                      f[f$variable == "huss"    & f$gcm == m,][[n, "file"]])
      }
    }

    else {
      ftas  <- f[f$variable == "tas"     & f$gcm == m,][[n, "file"]]
      frsds <- f[f$variable == "rsds"    & f$gcm == m,][[n, "file"]]
      fsfc  <- f[f$variable == "sfcwind" & f$gcm == m,][[n, "file"]]
      fhuss <- f[f$variable == "huss"    & f$gcm == m,][[n, "file"]]
    }

    print(paste("Processing temperature file:", ftas))

    if(bait) {
      hddcddCell <- calcStackHDDCDD(ftas,
                                    t_lim,
                                    countries,
                                    pop,
                                    hddcddFactor,
                                    bait,
                                    frsds = frsds,
                                    fsfc  = fsfc,
                                    fhuss = fhuss,
                                    wBAIT = wBAIT,
                                    params = params)}

    else {
      hddcddCell <- calcStackHDDCDD(ftas,
                                    t_lim,
                                    countries,
                                    pop,
                                    hddcddFactor,
                                    bait)}
    return(hddcddCell)
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



  # PARAMETERS------------------------------------------------------------------

  terra::setGDALconfig(c("BIGTIFF = YES"))


  # threshold temperature for heating and cooling [C]
  # NOTE: Staffel gives global average of T_heat = 14, T_cool = 20
  # t_lim <- list("HDD" = seq(12, 18), "CDD" = seq(20, 26))
  t_lim <- list("HDD" = seq(14), "CDD" = seq(22))

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

  # BAIT parameter names
  parNames <- c("a_rsds", "b_rsds", "a_sfcwind", "b_sfcwind", "a_huss", "b_huss")



  # READ-IN DATA----------------------------------------------------------------

  # list of files that are processed
  files <- toolGetMapping(mappingFile, type = "sectoral", where = "mappingfolder") %>%
    filter(variable != "")

  # cells -> country mask
  fCM <- file.path(files[files$variable == "CountryMask", "file"])
  countries <- readSource("ISIMIPbuildings", subtype = fCM, convert = FALSE)



  # PROCESS DATA----------------------------------------------------------------

  # calculate HDD/CDD-factors
  hddcddFactor <- calcHDDCDDFactors(tlow=-100.15, tup=74.85, t_lim, tamb_std, tlim_std)

  ssp   <- unique(files$ssp[files$variable == "tas"])
  rcp   <- unique(files$rcp[files$variable == "tas"])
  model <- unique(files$gcm[files$variable == "tas"])

  # read population data
  fpop <- files %>% filter(variable == "pop")
  pop  <- readSource("ISIMIPbuildings", subtype = fpop$file,
                     convert = FALSE)

  if (bait) {
    baitPars <- calcOutput("BAITpars", aggregate = FALSE, model = model)
    names(baitPars) <- parNames
  }

  hddcdd <- do.call( # file iteration
    "rbind",
    lapply(
      seq(nrow(filter(files, files$variable == "tas"))),
      function(n) {
        split <- files[files$variable == "tas" & files$gcm == model,][[n, "split"]]

        # split large raster files to save memory / speed up processing
        if (as.logical(split)) {
          hddcddCell <- do.call(
            "rbind",
            lapply(list("A", "B"), function(suffix) {
              tmp <- makeCalculations(f = files,
                                      m = model,
                                      n = n,
                                      t_lim = t_lim,
                                      countries = countries,
                                      pop = pop,
                                      hddcddFactor = hddcddFactor,
                                      bait = bait,
                                      wBAIT = wBAIT,
                                      params = baitPars,
                                      suffix = suffix)
              return(tmp)
            }
            )
          )
        }

        # smaller raster files
        else {
          hddcddCell <- makeCalculations(f = files,
                                         m = model,
                                         n = n,
                                         t_lim = t_lim,
                                         countries = countries,
                                         pop = pop,
                                         hddcddFactor = hddcddFactor,
                                         bait = bait,
                                         wBAIT = wBAIT,
                                         params = baitPars)
        }

        hddcddCell <- hddcddCell %>%
          mutate("model" = model,
                 "ssp" = ssp,
                 "rcp" = rcp)

        return(hddcddCell)
      }
    )
  )

  rownames(hddcdd) <- c()



  # OUTPUT----------------------------------------------------------------------

  # might be necessary for magpie output
  data <- hddcdd %>%
    mutate(period = gsub("y", "", .data[["period"]])) %>%
    unite("variable", .data[["variable"]], .data[["tlim"]], sep = "_") %>%
    unite("scenario", .data[["ssp"]], .data[["rcp"]], sep = "_")

  # prepare output
  data <- data %>%
    as.quitte() %>%
    as.magpie() %>%
    toolCountryFill()


  return(list(x = data))
}



