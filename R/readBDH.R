#' Read German heating system data from BDH
#'
#' @returns MagPIE object with yearly sales of heating systems in Germany
#'
#' @author Robin Hasse
#'
#' @importFrom utils read.csv
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate select
#' @importFrom magclass as.magpie collapseDim
#' @export

readBDH <- function() {
  read.csv("Marktentwicklung.csv", check.names = FALSE) %>%
    pivot_longer(c(-"period", -"gesamt"), names_to = "technology") %>%
    mutate(value = .data[["gesamt"]] * .data[["value"]]) %>%
    select(-"gesamt") %>%
    as.magpie(temporal = "period", datacol = "value") %>%
    collapseDim(1)
}
