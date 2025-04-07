calcUvalueEDGEB <- function(endOfHistory = 2025) {

  # READ-IN DATA ---------------------------------------------------------------

  ## U-value data ====

  library(stringr)   # str_to_lower
  library(MASS)      # boxcox
  library(dplyr)     # recode
  library(tibble)    # deframe

  hotmapsData <- readSource("Hotmaps") %>%
    as.quitte(na.rm = TRUE)

  messageData <- readSource("MessageIXBuildings", subtype = "uvalue") %>%
    as.quitte(na.rm = TRUE)

  euBuildObsData <- readSource("EUBuildingsDB", "BuildingShellPerformance") %>%
    as.quitte(na.rm = TRUE)


  # Floor space data ====

  messageFloorspace <- readSource("MessageIXBuildings", subtype = "floorByCohort") %>%
    as.quitte(na.rm = TRUE)

  odysseeData <- rbind(readSource("Odyssee", subtype = "households") %>%
                         as.quitte(na.rm = TRUE),
                       readSource("Odyssee", subtype = "services") %>%
                         as.quitte(na.rm = TRUE))


  # Heating / Cooling Degree-Days
  hddcdd <- readSource("HDDCDDtemp", convert = FALSE) %>%
    as.quitte()

  # Population
  pop <- calcOutput("Population", aggregate = FALSE) %>%
    as.quitte()

  # GDP per capita
  gdppop <- calcOutput("GDPpc",
                       scenario = "SSP2",
                       average2020 = FALSE,
                       unit = "constant 2005 Int$PPP",
                       aggregate = FALSE,
                       years = 1960:endOfHistory) %>%
    as.quitte()

  # TODO: use built-in weights of readHDDCDD



  # PARAMETERS -----------------------------------------------------------------

  # Minimum assumed U-Value
  uvalueMin <- 0.2



  ## Relevant parameters ====

  uvalueParameters <- c("BUILDING|Construction features (U-values)|Floor",
                        "BUILDING|Construction features (U-values)|Walls",
                        "BUILDING|Construction features (U-values)|Roof",
                        "BUILDING|Construction features (U-values)|Windows")

  buildingsToIgnore <- c("Residential sector|Total|Total",
                         "Residential sector|Residential sector|Total",
                         "Service sector|Total|Total")

  floorParameters <- c("BUILDING|Area|Constructed area")


  ## Mappings ====

  buildingTypeMap <- toolGetMapping("buildingTypeMapping_Hotmaps-BEAM2.csv",
                                   type = "sectoral",
                                   where = "mredgebuildings")

  buildingCharacteristics <- toolGetMapping("buildingCharacteristics_BEAM2.csv",
                                            type = "sectoral",
                                            where = "mredgebuildings")

  vintageMapMessage <- toolGetMapping("vintageMapping_MessageIX.csv",
                                      type = "sectoral",
                                      where = "mredgebuildings")

  varMappingEUBuildObs <- toolGetMapping("variableMapping_EUBuildingsObservatory.csv",
                                         type = "sectoral",
                                         where = "mredgebuildings")

  variableMap <- list("Cellar ceiling"       = "floor",
                      "Roof / upper ceiling" = "roof",
                      "Exterior Walls"       = "walls",
                      "Windows"              = "windows")

  variableMapOdyssee <- list("surter_m2" = "floorCom",        # total floor area of services in m2
                             "nbrlog_1"  = "nDwellings",      # number of residential dwellings
                             "surlog_m2" = "residentialAvg")  # average floor space of single dwelling in m2


  regionmapMessage <- toolGetMapping("regionmappingMessageIX.csv",
                               type = "regional", where = "mredgebuildings")



  # PROCESS DATA ---------------------------------------------------------------


  ## Process U-Value Data ====


  ### Hotmaps Data =====

  # U-Values per building vintage and building component (e.g. floor, windows, ...)
  uvaluesHotmapsData <- hotmapsData %>%
    filter(.data$variable %in% uvalueParameters,
           !.data$building %in% buildingsToIgnore) %>%
    mutate(variable = str_to_lower(word(variable, -1, sep = "\\|"))) %>%
    select(-"unit", -"model", -"scenario")

  # Prepare building components
  buildingCharacteristics <- buildingCharacteristics %>%
    mutate(variable = recode(.data$variable, !!!variableMap)) %>%
    rename("weight" = "value") %>%
    select(-"unit", -"source", -"floorspace", -"component_area")

  # Aggregate U-Values per building bage with weighted mean across component areas
  uvaluesAggComp <- uvaluesHotmapsData %>%
    left_join(buildingTypeMap, by = c("building" = "typHotmaps")) %>%
    left_join(buildingCharacteristics, by = c("variable", "typBEAM" = "building")) %>%
    group_by(across(all_of(c("region", "period", "bage", "building")))) %>%
    reframe(value = sum(.data$weight * .data$value)) %>%
    ungroup()

  ### SEPARATE HERE IN RESOLUTION

  # Floor space as aggregation weights
  floorspace <- hotmapsData %>%
    filter(.data$variable %in% floorParameters,
           !.data$building %in% buildingsToIgnore) %>%
    rename("floorspace" = "value") %>%
    select(-"scenario", -"model", -"unit", -"variable")

  # Aggregate u-values across building bages and types
  uvaluesHotmaps <- uvaluesAggComp %>%
    filter(!is.na(.data$value)) %>%
    left_join(floorspace, by = c("region", "period", "bage", "building")) %>%
    group_by(across(all_of(c("region", "period")))) %>%
    reframe(uvalue = sum(.data$value * .data$floorspace) / sum(.data$floorspace))




  ### EU Building Observatory ====

  # Floor space weights
  floorspaceOdyssee <- odysseeData %>%
    # rename variables
    filter(.data$variable %in% names(variableMapOdyssee)) %>%
    mutate(variable = recode(.data$variable, !!!variableMapOdyssee)) %>%
    select("region", "period", "variable", "value") %>%

    # calculate total residential floor space
    pivot_wider(names_from = "variable", values_from = "value") %>%
    mutate(floorRes = .data$nDwellings * .data$residentialAvg) %>%
    select(-"nDwellings", -"residentialAvg") %>%

    # calculate regional residential share
    group_by(across(all_of(c("region", "period")))) %>%
    mutate(floorShareRes = .data$floorRes / (.data$floorRes + .data$floorCom)) %>%
    ungroup() %>%

    # calculate global average shares per period
    group_by(across(all_of("period"))) %>%
    mutate(floorShareResGlo = mean(.data$floorShareRes, na.rm = TRUE)) %>%
    ungroup() %>%

    # allocate regional shares; use global average where NA
    mutate(floorShareRes = ifelse(is.na(.data$floorShareRes),
                                  .data$floorShareResGlo,
                                  .data$floorShareRes)) %>%
    select("region", "period", "floorShareRes")


  # Determine scaling factor component_area / floor space in (m2/m2)
  scalingFactor <- buildingCharacteristics %>%
    mutate(building = recode(.data$building,
                             !!!deframe(buildingTypeMap %>%
                                          select("typBEAM", "typHotmaps")))) %>%
    # join floor space
    right_join(floorspace,
               by = "building",
               relationship = "many-to-many") %>%
    filter(!is.na(.data$weight)) %>%

    # make variable declaration compliant with other data
    mutate(building = sub(" sector\\|.*", "", .data$building),
           building = ifelse(.data$building == "Residential", "residential", "commercial")) %>%

    # aggregate factors per building type and vintage
    group_by(across(all_of(c("region", "building", "variable")))) %>%
    reframe(weight = sum(.data$weight * .data$floorspace) / sum(.data$floorspace))


  # Calculate U-values per region and period
  uvaluesEUBuildObs <- euBuildObsData %>%
    filter(.data$variable %in% varMappingEUBuildObs$variable,
           .data$period != 2008) %>%

    left_join(varMappingEUBuildObs, by = "variable") %>%
    select("region", "period", "variable" = "variableNew", "building", "value") %>%

    left_join(scalingFactor, by = c("region", "variable", "building")) %>%

    group_by(across(all_of(c("region", "period", "building")))) %>%
    reframe(value = sum(.data$value * .data$weight)) %>%

    # make variable names compliant
    pivot_wider(names_from = "building", values_from = "value") %>%

    # aggregate u-values
    left_join(floorspaceOdyssee, by = c("region", "period")) %>%

    mutate(uvalue = .data$residential * .data$floorShareRes +
             .data$commercial * (1 - .data$floorShareRes)) %>%

    select("region", "period", "uvalue")



  ### MessageIX-Buildings / STURM ====

  # Pre-process floor space weights
  messageFloorspace <- messageFloorspace %>%
    filter(.data$scenario == "SSP2",
           .data$period <= endOfHistory) %>%
    select("region", "period", "floorspace" = "value", "variable" = "cohort")


  # Aggregate u-values across building vintages per region
  uvaluesMessage <- messageData %>%
    select("region", "variable", "value") %>%

    # aggregate to vintage types since u-values are not distinguished between single-/multi-family homes
    mutate(variable = sub("^[^_]+_", "", .data$variable)) %>%
    group_by(across(all_of(c("region", "variable")))) %>%
    reframe(uvalue = mean(.data$value)) %>%

    # match cohort types
    right_join(vintageMapMessage,
               by = c("variable" = "cohortCode"),
               relationship = "many-to-many") %>%

    # aggregate u-values by floor space
    right_join(messageFloorspace, by = c("region", "cohort" = "variable")) %>%
    group_by(across(all_of(c("region", "period")))) %>%
    reframe(uvalue = sum(.data$uvalue * .data$floorspace) / sum(.data$floorspace),
            floorspace = sum(.data$floorspace))



  ## Prepare data for regression ====

  # Aggregate u-values to MessageIX regions

  regionmapMessage <- regionmapMessage %>%
    select("region" = "CountryCode", "regionTarget" = "X11regions") %>%
    semi_join(uvaluesMessage, by = "region")

  uvaluesMessageAgg <- uvaluesMessage %>%
    mutate(unit = NA) %>%
    pivot_longer(cols = c("uvalue", "floorspace"), names_to = "variable", values_to = "value") %>%
    aggregate_map(mapping = regionmapMessage,
                  by = "region",
                  subset2agg = "uvalue",
                  weights = "floorspace",
                  forceAggregation = TRUE) %>%
    select(-"unit") %>%
    filter(variable == "uvalue") %>%
    pivot_wider(names_from = "variable", values_from = "value")


  # Filter correct historical limit temperatures
  # TODO: Needs to be made a lot better as soon as bette degree day data is there
  hddcdd <- hddcdd %>%
    filter(.data$variable %in% c("HDD_14", "CDD_20"),
           period != 201,
           !is.na(.data$value)) %>%
    mutate(variable = sub("_.*", "", .data$variable)) %>%
    filter(.data$period <= endOfHistory) %>%
    select("region", "period", "variable", "value") %>%
    group_by(across(all_of(c("region", "period", "variable")))) %>%
    reframe(value = mean(.data$value)) %>%
    pivot_wider(names_from = "variable", values_from = "value") %>%
    mutate(HDDCDD = .data$HDD + .data$CDD) %>%
    select(-"HDD", -"CDD") %>%
    pivot_longer(cols = c("HDDCDD"), names_to = "variable", values_to = "value")


  # Aggregate HDD/CDDS to MessageIX regions
  pop <- pop %>%
    filter(.data$variable == "pop_SSP2") %>%
    mutate(variable = sub("_SSP2", "", .data$variable)) %>%
    select("region", "period", "variable", "value")

  hddcddAgg <- hddcdd %>%
    rbind(pop) %>%
    mutate(unit = NA) %>%
    aggregate_map(mapping = regionmapMessage %>%
                    semi_join(hddcdd, by = "region"),
                  by = "region",
                  subset2agg = "HDDCDD",
                  weights = "pop",
                  forceAggregation = TRUE) %>%
    select(-"unit")

  # Full historical HDD/CDD data for all relevant regions
  hddcddHist <- hddcdd %>%
    filter(!.data$region %in% unique(regionmapMessage$regionTarget)) %>%
    rbind(hddcddAgg) %>%
    pivot_wider(names_from = "variable", values_from = "value")


  # Aggregate GDP/pop to MessageIX regions
  gdppop <- gdppop %>%
    mutate(variable = "gdppop") %>%
    select("region", "period", "variable", "value")

  gdppopAgg <- gdppop %>%
    rbind(pop) %>%
    mutate(unit = NA) %>%
    aggregate_map(mapping = regionmapMessage %>%
                    semi_join(gdppop, by = "region"),
                  by = "region",
                  subset2agg = "gdppop",
                  weights = "pop",
                  forceAggregation = TRUE) %>%
    select(-"unit")

  # Full historical GDP/pop data for all relevant regions
  gdppopHist <- gdppop %>%
    filter(!.data$region %in% unique(regionmapMessage$regionTarget)) %>%
    rbind(gdppopAgg) %>%
    pivot_wider(names_from = "variable", values_from = "value")


  # Gather full historical data on u-values, HDD/CDDs and GDP/pop
  dataHist <- uvaluesHotmaps %>%
    rbind(uvaluesEUBuildObs,
          uvaluesMessageAgg) %>%
    right_join(hddcddHist, by = c("region", "period"),
               relationship = "many-to-many") %>%
    left_join(gdppopHist, by = c("region", "period"))



  ## Make linear regression to fill historical data ====

  # Prune historical data to existing data points
  dataEstimate <- dataHist %>%
    filter(!is.na(.data$uvalue))


  ### Employ Box-Cox method to obtain optimal transformation parameter lambda ====

  # Shift U-values to ensure positivity
  dataEstimate$uvalueShifted <- dataEstimate$uvalue - uvalueMin

  # Perform Box-Cox transformation
  boxcoxEstimate <- lm("uvalueShifted ~ HDDCDD", data = dataEstimate) %>%
    boxcox(plotit = FALSE)

  # shift parameter
  lambda <- boxcoxEstimate$x[which.max(boxcoxEstimate$y)]

  # Apply Box-Cox transformation with the shift
  if (abs(lambda) < 1e-6) {
    dataEstimate$uvalueTrans <- log(dataEstimate$uvalueShifted)
  } else {
    dataEstimate$uvalueTrans <- (dataEstimate$uvalueShifted^lambda - 1) / lambda
  }


  # Fit a linear model on transformed data
  estimate <- lm("uvalueTrans ~ HDDCDD", data = dataEstimate)


  # Predict transformed values on full historical data set
  dataHist$predTrans = predict(estimate, newdata = dataHist)

  # Apply inverse Box-Cox transformation with shift
  dataHist <- dataHist %>%
    mutate(predUvalue = ifelse(abs(lambda) < 1e-6,
                               uvalueMin + exp(.data$predTrans),
                               uvalueMin + (lambda * .data$predTrans + 1)^(1 / lambda)))







}
