#' Calculate useful energy demand for space heating in buildings
#'
#' @author Robin Hasse
#'
#' @importFrom madrat readSource calcOutput toolGetMapping
#' @importFrom magclass mselect as.magpie collapseDim
#' @importFrom dplyr select right_join group_by across all_of summarise mutate
#'   cross_join
#' @importFrom quitte as.quitte
#' @importFrom brick getBrickMapping
#' @export

calcUEdemand <- function() {

  # map Hotmaps vintages
  vinMap <- toolGetMapping("vintageMapping_Hotmaps.csv", "sectoral",
                           "mredgebuildings") %>%
    select("vin", "vinHotmaps")

  # map Hotmaps building types
  typMap <- toolGetMapping("buildingTypeMapping_Hotmaps.csv", "sectoral",
                           "mredgebuildings") %>%
    select("typ", "typHotmaps")

  # Useful energy demand for space heating (kWh/yr/m2)
  ueDem <- readSource("Hotmaps") %>%
    mselect(variable = "ENERGY|Useful energy demand|Space heating") %>%
    toolCountryFillAvg(verbosity = 2) %>%
    as.quitte(na.rm = TRUE) %>%
    right_join(typMap, by = c(building = "typHotmaps")) %>%
    right_join(vinMap, by = c(bage = "vinHotmaps"),
               relationship = "many-to-many") %>%

    select("region", "typ", "vin", "value") %>%
    group_by(across(-all_of(c("value")))) %>%
    summarise(value = mean(.data[["value"]]), .groups = "drop")

  # add dimension: building shell
  relDem <- getBrickMapping("buildingShell.csv") %>%
    select("bs", relDem)
  ueDem <- ueDem %>%
    cross_join(relDem) %>%
    mutate(value = .data[["value"]] * .data[["relDem"]]) %>%
    select(-"relDem")

  # convert to magpie object
  ueDem <- ueDem %>%
    as.magpie(spatial = "region", datacol = "value") %>%
    collapseDim()

  # weight: FE demand
  feBuildings <- calcOutput("WeightFeBuildings", aggregate = FALSE) %>%
    mselect(period = "y2017", collapseNames = TRUE)

  return(list(x = ueDem,
              weight = feBuildings,
              unit = "kWh/yr/m2",
              min = 0,
              description = "Floor-space specific useful energy demand for space heating"))
}
