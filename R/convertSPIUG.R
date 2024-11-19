#' Convert heating systems stock and sales in Poland from SPIUG
#'
#' @param x MagPIE object with PIUG data
#' @returns MagPIE object with ISO 3 country codes
#'
#' @author Robin Hasse
#'
#' @importFrom madrat toolCountry2isocode toolCountryFill
#' @importFrom magclass add_dimension getItems getItems<-
#' @export

convertSPIUG <- function(x) {

  # add region dimension
  x <- add_dimension(x, 1, "region", "Poland")
  getItems(x, 1) <- toolCountry2isocode(getItems(x, 1))

  # fill missing countries with NA
  x <- toolCountryFill(x, verbosity = 2)

  return(x)
}
