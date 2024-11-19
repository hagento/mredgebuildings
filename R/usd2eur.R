#' Exchange rate between USD and EUR
#'
#' @param year integer, reference year
#' @returns MER EUR in USD in given year
#'
#' @author Robin Hasse
#'
#' @importFrom dplyr %>% .data filter

usd2eur <- function(year = 2020) {
  GDPuc:::wb_wdi %>% # nolint: undesirable_operator_linter.
    filter(.data[["year"]] == !!year, .data[["iso3c"]] == "DEU") %>%
    getElement("MER (LCU per US$)") %>%
    mean()
}
