#' Share of energy carriers or end uses in historic buildings demand
#'
#' Shares are calculated based on energy demands from the Odyssee data base for
#' EU member states,
#'
#' @author Robin Hasse, Antoine Levesque
#'
#' @param subtype Character, dimension names, level of shares.
#' @returns MAgPIE object with historic shares
#'
#' @importFrom madrat readSource toolGetMapping
#' @importFrom dplyr %>% filter mutate select group_by across all_of any_of
#'   .data
#' @importFrom quitte as.quitte removeColNa replace_column getRegs
#' @export

calcSharesBuildingDemand <- function(subtype = c("enduse", "carrier")) {

  # zero shares are filled with eps
  eps <- 0.005



  # READ DATA ------------------------------------------------------------------

  daioglou <- do.call("rbind", lapply(
    c("cooking", "space_heating", "water_heating", "lighting"),
    function(enduse) {
      readSource("Daioglou", paste0("share_", enduse)) %>%
        as.quitte() %>%
        filter(!is.na(.data[["value"]]),
               .data[["variable"]] != "other") %>%
        mutate(enduse = enduse)
    }
  ))
  odyssee <- calcOutput("ShareOdyssee", "enduse_carrier", aggregate = FALSE) %>%
    as.quitte() %>%
    filter(!is.na(.data[["value"]]))
  urbanshare <- calcOutput("UrbanPast", aggregate = FALSE)
  etpEurRegions <- toolGetMapping("regionmappingEDGE.csv",
                                  type = "regional", where = "mredgebuildings")
  regionalCorrespondence <- toolGetMapping("regionalCorrespondenceEDGECarrier.csv",
                                           type = "sectoral",
                                           where = "mredgebuildings")  %>%
    gather("enduse", "fillRegion", -"region") %>%
    rename(etpRegion = "region")



  # CALCULATE SHARES -----------------------------------------------------------

  # map carriers and normalise
  mapCarriers <- list(
    petrol  = c("kerosine", "lpg", "fuel oil"),
    biotrad = c("traditional fuelwood", "candles"),
    biomod  = c("improved fuelwood", "biogas"),
    natgas  = c("natural gas"),
    elec    = c("electricity"),
    coal    = c("coal"),
    heat    = c("remote heating")
  )
  mapCarriers <- do.call("rbind", lapply(names(mapCarriers), function(carrier) {
    data.frame(variable = mapCarriers[[carrier]], carrier = carrier)
  }))
  daioglou <- daioglou %>%
    replace_column(mapCarriers, "variable", "carrier") %>%
    rename(carrier = "variable") %>%
    group_by(across(all_of(
      c("region", "period", "quintile", "demographic", "carrier", "enduse")
    ))) %>%
    summarise(value = sum(.data[["value"]]), .groups = "drop") %>%
    ungroup() %>%
    mutate(value = ifelse(.data[["value"]] == 0, eps, .data[["value"]])) %>%
    group_by(across(all_of(
      c("region", "period", "quintile", "demographic", "enduse")
    ))) %>%
    mutate(value = .data[["value"]] / sum(.data[["value"]]))

  # aggregate demographics (population weighted) and drop quintiles
  urbanshare <- urbanshare %>%
    as.quitte() %>%
    select("region", "period", urbanshare = "value")
  daioglou <- daioglou %>%
    group_by(across(all_of(c("region", "period", "enduse")))) %>%
    filter(!any(.data[["demographic"]] == "Total"),
           .data[["quintile"]] == "0") %>%
    ungroup() %>%
    left_join(urbanshare, by = c("region", "period")) %>%
    group_by(across(all_of(c("region", "period", "enduse", "carrier", "quintile")))) %>%
    summarise(
      demographic = "Total",
      value = if (all(c("Rural", "Urban") %in% .data[["demographic"]])) {
        .data[["value"]][.data[["demographic"]] == "Urban"] *
          unique(.data[["urbanshare"]]) +
          .data[["value"]][.data[["demographic"]] == "Rural"] *
            (1 - unique(.data[["urbanshare"]]))
      } else {
        mean(.data[["value"]])
      },
      .groups = "drop"
    ) %>%
    rbind(daioglou) %>%
    filter(.data[["demographic"]] == "Total",
           .data[["quintile"]] == "0") %>%
    select(-"demographic", -"quintile")



  # AVERAGE SHARES IN WORLD REGIONS --------------------------------------------

  # average shares in world regions to fill missing countries
  daioglouAvg <- daioglou %>%
    replace_column(etpEurRegions, region = "CountryCode", "RegionCodeEUR_ETP") %>%
    group_by(across(all_of(c("region", "enduse", "carrier")))) %>%
    summarise(value = mean(.data[["value"]]), .groups = "drop") %>%
    ungroup() %>%
    mutate(value = ifelse(.data[["value"]] == 0, eps, .data[["value"]])) %>%
    group_by(across(all_of(c("region", "enduse")))) %>%
    mutate(value = .data[["value"]] / sum(.data[["value"]]))

  # space heating and water heating share assumed identical if one is missing
  daioglouAvg <- daioglouAvg %>%
    spread("enduse", "value") %>%
    mutate(space_heating = ifelse(is.na(.data[["space_heating"]]),
                                  .data[["water_heating"]],
                                  .data[["space_heating"]]),
           water_heating = ifelse(is.na(.data[["water_heating"]]),
                                  .data[["space_heating"]],
                                  .data[["water_heating"]])) %>%
    gather("enduse", "value", any_of(unique(daioglou$enduse))) %>%
    filter(!is.na(.data[["value"]]))

  # Complete lighting with electricity only
  daioglouAvg <- rbind(daioglouAvg, data.frame(region = "REP",
                                               enduse = "lighting",
                                               carrier = "elec",
                                               value = 1))

  # fill regions according to regional correspondence
  share <- etpEurRegions %>%
    select(region = "CountryCode", etpRegion = "RegionCodeEUR_ETP") %>%
    filter(!.data[["region"]] %in% getRegs(odyssee)) %>%
    left_join(regionalCorrespondence, by = "etpRegion") %>%
    left_join(daioglouAvg, by = c(fillRegion = "region", "enduse")) %>%
    rename(fillValue = "value") %>%
    left_join(daioglou, by = c("region", "carrier", "enduse"))

  share <- share # remove this line when continuing






}
