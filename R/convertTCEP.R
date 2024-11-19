#' Convert TCEP data base
#'
#' @param subtype character, type of data
#' @param x MAgPIE object with data from TCEP #nolint
#'
#' @returns magpie object
#'
#' @author Hagen Tockhorn, Robin Hasse
#'
#' @importFrom dplyr %>% select mutate
#' @importFrom rlang .data
#' @importFrom tidyr gather
#' @importFrom quitte as.quitte aggregate_map revalue.levels
#' @importFrom magclass as.magpie
#' @importFrom madrat toolGetMapping calcOutput
#' @export

convertTCEP <- function(x, subtype = "enduse") {

  if (subtype != "enduse") {
    stop("No conversion for this subtype: ", subtype)
  }

  # READ-IN DATA ---------------------------------------------------------------

  data <- x


  # Enduse Mapping
  tcepMap <- toolGetMapping("regionmappingTCEP.csv",
                            type = "regional",
                            where = "mredgebuildings")

  # Population
  # TODO: check if FE demand in buildings should be used as weights #nolint‚
  pop <- calcOutput("PopulationPast", aggregate = FALSE) %>%
    as.quitte()


  # PARAMETERS -----------------------------------------------------------------

  # enduse mapping
  enduseMap <- toolGetMapping(name = "enduseMap_TCEP.csv",
                              type = "sectoral",
                              where = "mredgebuildings") %>%
    pull("EDGE", "TCEP")


  # PROCESS DATA ---------------------------------------------------------------

  # remap enduse names
  data <- data %>%
    as.quitte() %>%
    select("region", "period", "variable", "unit", "value") %>%
    revalue.levels(variable = enduseMap)

  # correct FE enduse data for "appliances" in India
  data <- data %>%
    mutate(value = ifelse(.data[["region"]] == "India and other dev. Asia" &
                            .data[["period"]] == 2000 &
                            .data[["variable"]] == "appliances",
                          .data[["value"]] * 0.70,
                          .data[["value"]]))


  # Disaggregate to ISO Level
  data <- data %>%
    mutate(region = gsub("_", ".", .data[["region"]])) %>%
    aggregate_map(
      mapping = tcepMap[, c("CountryCode", "RegionCode")],
      by = c("region" = "RegionCode"),
      weights = pop %>%
        rename(weight_val_col = "value") %>%
        select("region", "period", "weight_val_col"),
      weight_item_col = "region",
      forceAggregation = TRUE
    ) %>%
    select("region", "period", "variable", "value") %>%
    rename(enduse = "variable")



  # OUTPUT ---------------------------------------------------------------------

  data <- data %>%
    as.magpie()

  return(data)

}
