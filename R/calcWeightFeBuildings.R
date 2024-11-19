#' calculate FE demand in buildings
#'
#' Simple estimate of total FE demand in different building types
#'
#' Should only be used as weights not actual demand. The disaggregation to
#' building types is too simplistic.
#'
#' @author Robin Hasse
#'
#' @importFrom dplyr %>%
#' @importFrom quitte inline.data.frame
#' @importFrom madrat readSource toolAggregate
#' @importFrom magclass mselect collapseDim getItems time_interpolate mbind
#'   setNames getItems
#' @export

calcWeightFeBuildings <- function() {

  # total FE demand in buildings subsectors
  fe <- readSource("IEA", "EnergyBalances") %>%
    mselect(FLOW = c("RESIDENT", "COMMPUB"),
            PRODUCT = "TOTAL") %>%
    collapseDim()

  # temporal interpolation
  fe[fe == 0] <- NA
  tfull <- getItems(fe, 2)
  fe <- do.call(mbind, lapply(getItems(fe, 1), function(r) {
    do.call(mbind, lapply(getItems(fe, 3.1), function(f) {
      t <- is.na(fe[r, , f])
      t <- getItems(t, 2)[t == 0]
      if (length(t) == 0) return(fe[r, tfull, f])
      time_interpolate(fe[r, t, f], tfull, extrapolation_type = "constant")
    }))
  }))
  fe <- toolCountryFillAvg(fe, verbosity = 2)

  # map building types to subsectors
  typMap <- inline.data.frame(
    "typ; FLOW;     weight",
    "SFH; RESIDENT; 0.5",
    "MFH; RESIDENT; 0.5",
    "Com; COMMPUB;  1.0"
  )
  weight <- typMap %>%
    group_by(across(all_of(c("typ", "weight")))) %>%
    reframe(region = getItems(fe, 1)) %>%
    ungroup() %>%
    as.magpie(spatial = "region", datacol = "weight") %>%
    collapseDim
  fe <- toolAggregate(fe,
                      rel = typMap, weight = weight,
                      from = "FLOW", to = "typ", dim = 3.1)
  getSets(fe) <- c("region", "period", "typ")

  return(list(x = fe,
              min = 0,
              unit = "EJ/yr",
              description = "Final energy demand of buildings (used as weights)"))
}
