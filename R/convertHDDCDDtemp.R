#' Convert temporary HDD and CDD data
#'
#' The temporary data already is at EDGE-B regional resolution. It is uniformly
#' disaggregated to country level here to comply with the madrat framework
#' before it is aggregated to the identical resolution in fullEDGEBUILDINGS.
#'
#' @author Robin Hasse
#'
#' @param x MAgPIE object with HDD and CDD data at EDGE-B resolution
#' @returns HDD CDD data at country resolution
#'
#' @importFrom madrat toolGetMapping toolAggregate
#' @export

convertHDDCDDtemp <- function(x) {
  mapping <- toolGetMapping(name = "regionmappingEDGE.csv",
                            type = "regional",
                            where = "mredgebuildings")

  toolAggregate(x, mapping,
                from = "RegionCodeEUR_ETP", to = "CountryCode")
}
