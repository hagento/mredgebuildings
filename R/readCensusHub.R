#' Read data from the European Census Hub
#'
#' The Census Hub of the European Statistical System provides census data from
#' 2011 on population and housing of EU member states. 2021 census data is
#' announced to be published there too.
#' - typeVintage: number of dwellings of different construction periods and
#' building types (single, two, multiple dwellings). Show data on: dwellings,
#' location: nations, Topics: Type of building, Period of construction
#' (select all) -> CSV (Separator: Comma)
#'
#' @source https://ec.europa.eu/CensusHub2/
#'
#' @param subtype census subset
#' @returns MAgPIE object with data
#'
#' @author Robin Hasse
#'
#' @importFrom madrat toolSubtypeSelect
#' @importFrom dplyr %>% select rename
#' @importFrom utils read.csv
#' @importFrom quitte as.quitte
#' @importFrom magclass as.magpie
#' @export

readCensusHub <- function(subtype) {

  # pick file
  files <- list(typeVintage = "csvoutput_HC53_2022_04_08_17_01.csv")

  data <- toolSubtypeSelect(subtype, files) %>%
    read.csv() %>%
    select(-"FLAGS", -"FOOTNOTES") %>%
    rename(region = "GEO",
           period = "TIME",
           value = "VALUE") %>%
    as.quitte() %>%
    as.magpie()

  return(data)
}
