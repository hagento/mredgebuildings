#' Convert heating systems sales in Italy from Assotermica
#'
#' @param x MagPIE object with Assotermica data
#' @returns MagPIE object with ISO 3 country codes
#'
#' @author Robin Hasse
#'
#' @importFrom madrat toolCountry2isocode toolCountryFill toolGetMapping
#' @importFrom magclass add_dimension getItems getItems<-
#' @export

convertAssotermica <- function(x) {

  # add region dimension
  x <- add_dimension(x, 1, "region", "Italy")
  getItems(x, 1) <- toolCountry2isocode(getItems(x, 1))

  # fill missing countries with NA
  x <- toolCountryFill(x, verbosity = 2)

  # translate technologies to English
  mapAssotermica <- toolGetMapping("technologyMapping_Assotermica.csv",
                                   type = "sectoral", where = "mredgebuildings")
  mapAssotermica <- stats::setNames(mapAssotermica[["variable"]],
                                    mapAssotermica[["variableITA"]])
  getItems(x, 3.1) <- unname(mapAssotermica[getItems(x, 3.1)])

  return(x)
}
