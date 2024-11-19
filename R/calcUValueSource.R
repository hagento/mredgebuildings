#' Passes U-values read from source
#'
#' @author Falk Benke, Robin Hasse
#'
#' @param subtype source of U-values
#' @export

calcUValueSource <- function(subtype = c("EUBuildingsObservatory", "ETSAP")) {

  subtype <- match.arg(subtype)

  x <- if (subtype == "EUBuildingsObservatory") {
    readSource("EUBuildingsObservatory", convert = FALSE)
  } else {
    readSource("ETSAP", convert = FALSE)
  }

  return(list(x = x,
              min = 0,
              isocountries = FALSE,
              unit = "W/m2K",
              description = paste("U-values as reported by", subtype)))
}
