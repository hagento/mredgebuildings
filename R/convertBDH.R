#' Convert German heating system data from BDH
#'
#' @param x MagPIE object with BDH data
#' @returns MagPIE object with ISO 3 country codes
#'
#' @author Robin Hasse
#'
#' @importFrom madrat toolCountry2isocode toolCountryFill toolGetMapping
#' @importFrom magclass add_dimension getItems getItems<-
#' @export

convertBDH <- function(x) {

  # add region dimension
  x <- add_dimension(x, 1, "region", "Germany")
  getItems(x, 1) <- toolCountry2isocode(getItems(x, 1))

  # fill missing countries with NA
  x <- toolCountryFill(x, verbosity = 2)

  # translate technologies to English
  mapBDH <- toolGetMapping("technologyMapping_BDH.csv", type = "sectoral",
                           where = "mredgebuildings")
  mapBDH <- stats::setNames(mapBDH[["technology"]], mapBDH[["technologyDEU"]])
  getItems(x, 3) <- unname(mapBDH[getItems(x, 3)])

  return(x)
}
