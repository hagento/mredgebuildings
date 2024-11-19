#' convert Hotmaps
#'
#' @author Robin Hasse
#'
#' @param x MagPIE object with raw building stock data
#' @returns MagPIE object with converted data
#'
#' @importFrom magclass getItems<- getItems
#' @importFrom madrat toolCountry2isocode toolCountryFill

convertHotmaps <- function(x) {

  data <- x

  # fix typos
  fixTypo <- function(m, reg, dim, from, to) {
    fix <- m[reg, , grep(from, getItems(x, dim), value = TRUE)]
    m <- m[, , grep(from, getItems(x, dim), value = TRUE), invert = TRUE]
    getItems(fix, dim) <- sub(from, to, getItems(fix, dim))
    m[reg, , getItems(fix, 3, full = TRUE)] <- fix
    return(m)
  }
  data <- data %>%
    fixTypo("cy", "bage", "Berfore", "Before") %>%
    fixTypo("pl", "variable", "\\(U-value\\)", "(U-values)") %>%
    fixTypo("ee", "building", "^\\|Appartment blocks", "Residential sector|Appartment blocks")

  # rename regions: ISO2 -> ISO3
  data <- data["eu27+UK", , invert = TRUE]
  getItems(data, 1) <- toolCountry2isocode(getItems(data, 1))

  # fill missing regions with NA
  data <- toolCountryFill(data, verbosity = 2)

  return(data)
}
