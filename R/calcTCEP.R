#' IEA's Tracking clean energy progress 2014
#'
#' @param subtype character, type of data
#' @author Robin Hasse
#'
#' @importFrom madrat readSource
#' @export

calcTCEP <- function(subtype = "enduse") {
  x <- switch(subtype,
    enduse = readSource("TCEP", "enduse", convert = FALSE),
    floorspace = readSource("TCEP", "floorspace", convert = FALSE)
  )
  return(list(x = x,
              min = 0,
              isocountries = FALSE,
              unit = "m2/cap",
              description = "Floor space per capita"))
}
