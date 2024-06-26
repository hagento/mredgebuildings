#' Converts Data from Primary Fuel and Useful Energy Data Base (PFUDB)
#'
#' PFU Data is disaggregated into different thermal and non-thermal Enduses
#' using FE Shares from Odyssee and IEA_ETP.
#'
#' This was adapted from EDGE function 'getPFUDB.R'.
#'
#' @returns magpie object
#'
#' @author Hagen Tockhorn, Robin Hasse
#'
#' @importFrom madrat readSource toolGetMapping readSource
#' @importFrom quitte aggregate_map as.quitte factor.data.frame revalue.levels
#' @importFrom dplyr %>% .data filter group_by summarise mutate ungroup select
#'   all_of any_of pull
#' @export

calcPFUDB <- function() {

  # FUNCTIONS ------------------------------------------------------------------

  # Sum to Carrier Level
  sumDf <- function(df, variables, newname) {
    carrierSum <- df %>%
      filter(.data[["carrier"]] %in% variables) %>%
      group_by(across(-any_of(c("value", "carrier")))) %>%
      summarise(value = sum(.data[["value"]], na.rm = TRUE), .groups = "drop") %>%
      mutate(carrier = newname) %>%
      ungroup() %>%
      select(all_of(colnames(df)))
    df %>%
      filter(!(.data[["carrier"]] %in% variables)) %>%
      rbind(carrierSum)
  }

  # replace NA's only if not all values are NA
  replaceNAwithZero <- function(x) {
    if(all(is.na(x))) {
      return(NA)
    } else {
      return(replace(x, is.na(x), 0))
    }
  }

  # convert enduse "appliances" to "refrigerators" with respective scaling factor
  addThermal <- function(df, mapping, fridgeShare, feOnly) {
    df <- df %>%
      filter(.data[["enduse"]] != "lighting") %>%
      left_join(regmappingEDGE %>%
                  select(-"RegionCodeEUR", -"RegionCodeEUR_ETP", -"X") %>%
                  rename(region = "CountryCode") %>%
                  left_join(fridgeShare, by = "RegionCode") %>%
                  select(-"RegionCode"),
                by = "region") %>%
      mutate(value = ifelse(.data[["enduse"]] != "appliances",
                            .data[["value"]],
                            .data[["value"]] * .data[["share"]]),
             enduse = ifelse(.data[["enduse"]] == "appliances",
                             "refrigerators",
                             as.character(.data[["enduse"]]))) %>%
      select(-"share")

    if (!feOnly) {
      df <- normalize(df, shareOf)
      return(df)
    }

    return(df)
  }



  # PARAMETERS -----------------------------------------------------------------

  # Rename Carrier Items
  carriersnames <- toolGetMapping("carrierMap_PFUDB.csv", "sectoral",
                                  "mredgebuildings") %>%
    pull("EDGE", "PFUDB")

  # Enduse-Carrier combinations which will be systematically excluded
  exclude <- toolGetMapping("excludeEnduseCarrier.csv", "sectoral",
                            "mredgebuildings")

  # fridge electricity shares (see calcShares)
  fridgeShare <- rbind(
    data.frame(RegionCode = "USA", share  = 0.12),
    data.frame(RegionCode = c("EUR", "OCD", "RUS", "JPN"), share = 0.17),
    data.frame(RegionCode = c("CHN", "IND", "NCD", "AFR", "MIE", "OAS"), share = 0.3))

  # lower temporal threshold of historical data
  periodBegin <- 1990



  # READ-IN DATA ---------------------------------------------------------------

  pfu <- readSource("PFUDB") %>%
    as.quitte()


  # Relevant toolDisaggregate Input-Data

  sharesEU <- calcOutput("Shares",
                         subtype = "enduse_thermal",
                         aggregate = TRUE,
                         regionmapping = "regionmappingEUshares.csv") %>%
    as.quitte()


  feOdyssee <- calcOutput("ShareOdyssee",
                          subtype = "enduse_carrier",
                          feOnly = TRUE,
                          aggregate = FALSE) %>%
    as.quitte()


  sharesOdyssee <- calcOutput("ShareOdyssee",
                          subtype = "enduse_carrier",
                          feOnly = FALSE,
                          aggregate = FALSE) %>%
    as.quitte()


  # FE IEA EEI data
  feIEAEEI <- calcOutput("IEA_EEI",
                         subtype = "buildings",
                         aggregate = FALSE) %>%
    as.quitte()



  # Mappings

  # EU Shares region mapping
  regmapping <- toolGetMapping(name  = "regionmappingEUshares.csv",
                               type  = "regional",
                               where = "mredgebuildings") %>%
    select("region", "regionAgg")

  # Enduse Mappings
  enduseMapping <- toolGetMapping(name  = "enduseMap_PFUDB.csv",
                                  type  = "sectoral",
                                  where = "mredgebuildings")

  # EDGE mapping
  regmappingEDGE <- toolGetMapping(name  = "regionmappingEDGE.csv",
                                   type  = "regional",
                                   where = "mredgebuildings")



  # PROCESS DATA ---------------------------------------------------------------

  pfu <- pfu %>%
    # Generalize Heat Carriers
    sumDf(c("Heat", "Geothermal", "Solar"), "heat") %>%

    # Map Carrier Names and Convert Units
    revalue.levels(carrier = carriersnames) %>%
    factor.data.frame() %>%
    mutate(value = replace_na(.data[["value"]], 0)) %>%
    select("region", "period", "unit", "carrier", "enduse", "value")  %>%
    semi_join(sharesEU, by = c("period"))



  ## Disaggregate into Thermal and Non-Thermal Part ====

  # Aggregate uses from PFU to appliances_light and keep the df as Non-thermal
  pfuNonTherm <- pfu %>%
    filter(.data[["enduse"]] != "Low-T heat") %>%
    aggregate_map(mapping = enduseMapping,
                  by = "enduse",
                  variable = "enduse",
                  forceAggregation = TRUE)


  ## Prepare toolDisaggregate Input ====

  # calculate enduse-carrier shares for IEA EEI data
  sharesIEAEEI <- feIEAEEI %>%
    group_by(across(all_of(c("region", "period")))) %>%
    mutate(value = replaceNAwithZero(.data[["value"]]),
           value = proportions(.data[["value"]])) %>%
    ungroup()

  sharesReplace <- sharesOdyssee %>%
    left_join(sharesIEAEEI, by = c("region", "period", "carrier", "enduse")) %>%
    mutate(value = ifelse(is.na(.data[["value.x"]]),
                          .data[["value.y"]],
                          .data[["value.x"]])) %>%
    select("region", "period", "carrier", "enduse", "value")

  # note: hard-coded the fridge share so that the FE data becomes compliant with
  # the remaining input (can be done more nicely)

  # correct fridgeShare
  sharesReplace <- sharesReplace %>%
    addThermal(regmappingEDGE, fridgeShare, feOnly = TRUE)

  feDisagg <- feOdyssee %>%
    left_join(feIEAEEI,
              by = c("region", "period", "carrier", "enduse")) %>%
    mutate(value = ifelse(is.na(.data[["value.x"]]),
                          .data[["value.y"]],
                          .data[["value.x"]])) %>%
    addThermal(regmappingEDGE, fridgeShare, feOnly = TRUE) %>%
    select("region", "period", "carrier", "enduse", "value")


  # Extract regions with existing disaggregated FE shares
  replaceRegs <- sharesReplace %>%
    filter(!is.na(.data[["value"]])) %>%
    pull("region") %>%
    unique()

  # re-aggregate Odyssee shares to carrier level
  sharesReplace <- sharesReplace %>%
    filter(.data[["region"]] %in% replaceRegs) %>%
    mutate(value = ifelse(is.na(.data[["value"]]),
                          0,
                          .data[["value"]])) %>%
    filter(.data[["enduse"]] != "lighting") %>%
    mutate(value = .data[["value"]] * ifelse(.data[["enduse"]] != "appliances",
                                             1,
                                             fridgeShare),
           enduse = ifelse(.data[["enduse"]] == "appliances",
                           "refrigerators",
                           as.character(.data[["enduse"]]))) %>%
    group_by(across(all_of(c("region", "period", "carrier")))) %>%
    mutate(value = proportions(.data[["value"]])) %>%
    ungroup()

  # exclude specific carrier/enduse combinations
  sharesReplace <- rbind(
    sharesReplace %>%
      semi_join(exclude, by = c("enduse", "carrier")) %>%
      mutate(value = 0),
    sharesReplace %>%
      anti_join(exclude, by = c("enduse", "carrier"))
  )


  sharesEU <- sharesEU %>%
    select("region", "period", "enduse", "value")


  # data excluded from disaggregation
  pfuThermExclude <- pfu %>%
    # exclude regions with existing disaggregated shares
    filter(.data[["enduse"]] == "Low-T heat",
           .data[["region"]] %in% replaceRegs) %>%
    select(-"enduse") %>%

    # disaggregate with existing shares
    left_join(sharesReplace %>%
                select("region", "period", "carrier", "enduse", "value") %>%
                rename(share = "value"),
              by = c("region", "period", "carrier"),
              relationship = "many-to-many") %>%
    mutate(value = .data[["value"]] * .data[["share"]],
           value = replace_na(.data[["value"]], 0)) %>%
    select(-"share")


  # Disaggregate Low-T Heat into different enduses
  pfuThermFE <- pfu %>%
    filter(.data[["enduse"]] == "Low-T heat",
           !(.data[["region"]] %in% replaceRegs),
           unit == "fe") %>%
    select(-"enduse") %>%
    toolDisaggregate(enduseShares  = sharesEU,
                     outliers      = c("IND", "CHN", "ZAF"),
                     exclude       = exclude,
                     dataDisagg    = feDisagg,
                     regionMapping = regmapping) %>%
    select(colnames(pfuNonTherm))

  # Use carrier-enduse distribution to apply on useful energy
  shares <- pfuThermFE %>%
    group_by(across(all_of(c("region", "period", "carrier")))) %>%
    mutate(share = proportions(.data[["value"]])) %>%
    select(-"value", -"unit")

  # disaggregate useful energy with calculated shares
  pfuTherm <- pfu %>%
    filter(.data[["enduse"]] == "Low-T heat",
           .data[["unit"]]   == "ue",
           !(.data[["region"]] %in% replaceRegs)) %>%
    select(-"enduse") %>%
    left_join(shares, by = c("region", "period", "carrier")) %>%
    mutate(value = .data[["value"]] * .data[["share"]]) %>%
    select(-"share") %>%
    rbind(pfuThermFE)


  # Join Non-Thermal and Thermal Part
  pfuRes <- rbind(pfuTherm, pfuThermExclude, pfuNonTherm)


  # Include "refrigerators" in "appliances_light"
  pfuRes <- rbind(
    pfuRes %>%
      filter(!(.data[["enduse"]] %in% enduseMapping[4:5, "enduse"])),
    pfuRes %>%
      aggregate_map(mapping = enduseMapping[4:5, ],
                    by = "enduse",
                    variable = "enduse",
                    forceAggregation = TRUE))

  # Select data above lower history boundary
  pfuRes <- pfuRes %>%
    filter(.data[["period"]] >= periodBegin)



  # OUTPUT ---------------------------------------------------------------------

  result <- pfuRes %>%
    as.quitte() %>%
    as.magpie()


  data <- list(
    x = result,
    weight = NULL,
    unit = "EJ",
    min = 0,
    description = "Primary Fuel and Useful Energy Data Base"
  )

  return(data)

}
