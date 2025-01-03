#' Historic and projected population
#'
#' Population differentiated by location (urban/rural) and building type
#'
#' @author Robin Hasse
#'
#' @importFrom magclass add_dimension mbind getItems getItems<- getSets<-
#'   time_interpolate getSets
#' @importFrom madrat calcOutput
#' @importFrom utils tail
#' @export

calcPopulationBuildings <- function() {

  pop <- calcOutput("Population", aggregate = FALSE)
  pop[pop < 0] <- 0
  getSets(pop) <- c("region", "period", "scenario")
  getItems(pop, "scenario") <- sub("pop_", "", getItems(pop, "scenario"))

  urbanShare <- calcOutput("Urban", aggregate = FALSE)
  getSets(urbanShare) <- c("region", "period", "scenario")
  getItems(urbanShare, "scenario") <- sub("urb_", "", getItems(urbanShare, "scenario"))

  # assume building type shares according to
  # TODO: improve split by building type # nolint: todo_comment_linter
  typShare <- calcOutput("BuildingStock", aggregate = FALSE)[, , "dwellings"]
  typShare <- dimSums(typShare, setdiff(tail(getSets(typShare, -2)), "typ"), na.rm = TRUE)
  typShare <- typShare / dimSums(typShare)
  typShare <- toolCountryFillAvg(typShare, verbosity = 2)

  # common periods
  t <- union(getItems(pop, 2), getItems(urbanShare, 2))
  pop <- time_interpolate(pop, t)
  urbanShare <- time_interpolate(urbanShare, t)
  typShare <- time_interpolate(typShare, t, extrapolation_type = "constant")

  # split population in rural/urban
  pop <- mbind(
    add_dimension(pop * urbanShare,       3.2, "loc", "urban"),
    add_dimension(pop * (1 - urbanShare), 3.2, "loc", "rural")
  )

  # split population in SFH/MFH
  pop <- pop * typShare

  return(list(x = pop,
              weight = NULL,
              unit = "million",
              description = "Rural and urban population",
              min = 0))
}
