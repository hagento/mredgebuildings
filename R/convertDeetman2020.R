#' Convert Weibull lifetime distribution parameters for buildings
#'
#' @param x MAgPIE object with data from Deetman et al. 2020
#' @param subtype character, building subsector (either 'residential' or
#'   'commercial')
#'
#' @author Robin Hasse
#'
#' @importFrom madrat toolGetMapping toolAggregate
#' @export

convertDeetman2020 <- function(x, subtype) {

  regionmapping <- toolGetMapping("regionmappingDeetman2020.csv",
                                  type = "regional", where = "mredgebuildings")

  switch(subtype,
    residential = toolAggregate(x, regionmapping),
    commercial = toolCountryFillAvg(x, verbosity = 2, no_remove_warning = "GLO")
  )
}
