#' Calculate Reference data for input matching
#'
#' @author Robin Hasse
#'
#' @param subtype character, matching reference
#'
#' @importFrom magclass mbind as.magpie collapseDim
#' @importFrom madrat readSource toolCountryFill toolGetMapping
#' @importFrom quitte as.quitte interpolate_missing_periods
#' @importFrom dplyr group_by filter mutate .data
#' @export
#'
calcMatchingReference <- function(subtype) {

  periods <- 2000:2020

  weight <- NULL
  minVal <- 0
  maxVal <- NULL
  description <- NULL

  # reference mapping
  refMap <- toolGetMapping(paste0("refMap_", subtype, ".csv"),
                           "sectoral", "mredgebuildings",
                           returnPathOnly = TRUE) %>%
    read.csv(comment.char = "#")



  # READ DATA ------------------------------------------------------------------

  if (grepl("Odyssee", subtype)) {
    readOdysseeData <- function(subsector) {
      readSource("Odyssee", subsector) %>%
        as.quitte(na.rm = TRUE) %>%
        select("region", "period", "variable", "value") %>%
        mutate(variable = sub("_.*$", "", .data[["variable"]]))
    }
    households <- readOdysseeData("households")
    services <- readOdysseeData("services")
  } else if (grepl("EUBDB", subtype)) {
    eubdbCateories <- c(
      "BuildingStockCharacteristics",
      "TechnicalBuildingSystems"
    )
    # labels cannot be longer than 63 characters in gams -> shorten
    eubdb <- do.call(rbind, lapply(eubdbCateories, function(c) {
      readSource("EUBuildingsDB", c) %>%
        as.quitte(na.rm = TRUE)
    })) %>%
      mutate(variable = sub("_.*$", "", .data[["variable"]]),
             variable = sub("Number", "No", .data[["variable"]]),
             variable = sub("non.residential", "com", .data[["variable"]]),
             variable = sub("dwelling", "dwel", .data[["variable"]])) %>%
      select("region", "period", "variable", "value")

  } else if (grepl("mredgebuildings", subtype)) {
    stockHist <- calcOutput("BuildingStock", aggregate = FALSE) %>%
      as.quitte(na.rm = TRUE) %>%
      filter(.data[["variable"]] == "floor") %>%
      select(-"model", -"scenario", -"variable", -"unit")
  }

  aggStock <- function(stockHist, dim) {
    data <- stockHist %>%
      group_by(across(all_of(c("region", "period", dim)))) %>%
      summarise(value = sum(.data[["value"]]) / 1E6,
                .groups = "drop")
    colnames(data)[3] <- "variable"
    return(as.quitte(data))
  }



  # PREPARE DATA ---------------------------------------------------------------

  # nolint start: todo_comment_linter.

  switch(subtype,

    ## mredgebuildings_buildingType ====

    mredgebuildings_buildingType = {

      data <- aggStock(stockHist, "buildingType")

      unit <- "million m2"
      description <- "stock of residential buildings"

    },

    ## mredgebuildings_location ====

    mredgebuildings_location = {

      data <- aggStock(stockHist, "location")

      unit <- "million m2"
      description <- "stock of residential buildings"

    },

    ## mredgebuildings_vintage ====

    mredgebuildings_vintage = {

      data <- aggStock(stockHist, "vintage")

      unit <- "million m2"
      description <- "stock of residential buildings"

    },

    ## mredgebuildings_heating ====

    mredgebuildings_heating = {

      data <- aggStock(stockHist, "heating")

      unit <- "million m2"
      description <- "stock of residential buildings"

    },

    ## IDEES_heatingShare ====

    IDEES_heatingShare = {

      data <- c("Residential", "Tertiary") %>%
        lapply(readSource, type = "JRC_IDEES") %>%
        do.call(what = mbind) %>%
        as.quitte(na.rm = TRUE) %>%
        select(-"scenario", -"model", -"unit")

      data <- refMap %>%
        select("variable", ".variableOriginal", ".variableTotal") %>%
        unique() %>%
        left_join(data,
                  by = c(.variableOriginal = "variable")) %>%
        group_by(across(-all_of(c(".variableOriginal", "value")))) %>%
        summarise(value = sum(.data[["value"]]), .groups = "drop") %>%
        left_join(data,
                  by = c(.variableTotal = "variable", "region", "period")) %>%
        mutate(value = .data[["value.x"]] / .data[["value.y"]],
               unit = "1") %>%
        select(-"value.x", -"value.y", -".variableTotal") %>%
        as.quitte()

      unit <- "1"
      description <- "share of heating systems in the stock"

    },

    ## Odyssee_stock ====

    Odyssee_stock = {

      vars <- list(
        households = unique(refMap[refMap[[".dataSet"]] == "households", "variable"]),
        services = unique(refMap[refMap[[".dataSet"]] == "services", "variable"])
      )
      households <- households %>%
        filter(.data[["variable"]] %in% vars$households) %>%
        mutate(value = .data[["value"]] / 1E6,
               unit = "million")
      services <- services %>%
        filter(.data[["variable"]] %in% c(vars$services)) %>%
        mutate(value = .data[["value"]] / 1E6,  # m2 -> million m2
               unit = "million m2")

      data <- rbind(households, services) %>%
        as.quitte()

      unit <- "million m2 or million dwellings"
      description <- "stock of dwellings and service floor space"

    },


    ## Odyssee_dwelSize ====

    Odyssee_dwelSize = {

      data <- refMap %>%
        select("variable") %>%
        unique() %>%
        left_join(households, by = "variable")

      # TODO: weight should be floor area -> circle dependency

      unit <- "m2/dwel"
      description <- "average dwelling size of the stock and new construction"

    },


    ## Odyssee_construction ====

    Odyssee_construction = {

      data <- refMap %>%
        select("variable") %>%
        unique() %>%
        left_join(households, by = "variable") %>%
        mutate(value = .data[["value"]] / 1E6,    # dwel/yr -> million dwel/yr
               unit = "million dwel/yr") %>%
        as.quitte()

      unit <- "million dwel/yr"
      description <- "flow of newly constructed buildings"

    },


    ## Odyssee_constructionFloor ====

    Odyssee_constructionFloor = {

      data <- refMap %>%
        select("variable") %>%
        unique() %>%
        separate("variable", c("dwel", "size"), sep = "_") %>%
        left_join(households, by = c(dwel = "variable")) %>%
        left_join(households, by = c(size = "variable", "region", "period")) %>%
        mutate(value = .data[["value.x"]] * .data[["value.y"]] / 1E6, # m2/yr -> million m2/yr
               unit = "million m2/yr") %>%
        select(-"value.x", -"value.y") %>%
        unite("variable", "dwel", "size", sep = "_") %>%
        as.quitte() %>%
        interpolate_missing_periods(periods, expand.values = TRUE) # TODO: remove extrapolation

      unit <- "million m2/yr"
      description <- "flow of newly constructed buildings"

    },



    ## Odyssee_heatingShare ====

    Odyssee_heatingShare = {
      data <- rbind(households, services)

      data <- refMap %>%
        select("variable") %>%
        unique() %>%
        separate("variable", c("heating", "total"), sep = "_") %>%
        left_join(households, by = c(heating = "variable")) %>%
        left_join(households, by = c(total = "variable", "region", "period")) %>%
        mutate(value = .data[["value.x"]] / .data[["value.y"]],
               unit = "1") %>%
        select(-"value.x", -"value.y") %>%
        unite("variable", "heating", "total", sep = "_") %>%
        as.quitte() %>%
        interpolate_missing_periods(periods) # TODO: remove extrapolation

      # TODO: weights

      unit <- "1"
      description <- "share of heating systems in the stock"

    },


    ## EUBDB_vintage ====

    EUBDB_vintage = {

      data <- refMap %>%
        select("variable") %>%
        unique() %>%
        left_join(eubdb, by = "variable") %>%
        mutate(unit = "1") %>%
        as.quitte()

      unit <- "1"
      description <- "share of dwellings in vintage cohort"

    },


    ## EUBDB_stock ====

    EUBDB_stock = {

      data <- refMap %>%
        select("variable", "qty") %>%
        unique() %>%
        left_join(eubdb, by = "variable") %>%
        mutate(value = .data[["value"]] / 1E6,    # qty/yr-> million qty/yr
               unit = c(dwel = "million dwel", area = "million m2")[.data[["qty"]]]) %>%
        select(-"qty") %>%
        as.quitte()

      unit <- "million dwel/yr"
      description <- "stock of dwellings and service floor space"
    },


    ## EuropeanCommissionRenovation ====

    EuropeanCommissionRenovation = {

      # The data represents the average between 2014 and 2017. But we assume
      # this value for all historic periods temporarily
      data <- readSource("EuropeanCommissionRenovation") %>%
        as.quitte(na.rm = TRUE) %>%
        unite("variable", "variable", "renovation", "subsector", sep = "|") %>%
        right_join(refMap %>%
                     select("variable") %>%
                     unique(),
                   by = "variable") %>%
        select(-"model", -"scenario") %>%
        group_by(across(all_of(c("region", "variable", "unit", "value")))) %>%
        complete(period = periods) %>%
        as.quitte()

      unit <- "1/yr"
      description <- "stock of dwellings and service floor space"
    },

    stop("The subtype '", subtype, "' is an invalid matching reference.")
  )

  # nolint end: todo_comment_linter.

  data <- data %>%
    filter(.data[["period"]] %in% periods) %>%
    as.magpie() %>%
    collapseDim(keepdim = c("region", "period", "variable")) %>%
    toolCountryFill(NA, verbosity = 2)


  return(list(x = data,
              weight = weight,
              unit = unit,
              min = minVal,
              max = maxVal,
              description = description))
}
