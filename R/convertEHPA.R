#' Correct heat pump sales data from EHPA
#'
#' @param x MagPIE object with EHPA data
#' @returns MagPIE object with ISO 3 country codes
#'
#' @author Robin Hasse
#'
#' @importFrom madrat toolCountry2isocode toolCountryFill
#' @export

convertEHPA <- function(x) {

  # common set names
  getSets(x)[1:2] <- c("region", "period")

  # rename regions: ISO2 -> ISO3
  x <- x["Total", , , invert = TRUE]
  getItems(x, 1) <- toolCountry2isocode(getItems(x, 1))

  # fill missing countries with NA
  x <- toolCountryFill(x, verbosity = 2)

  return(x)
}
