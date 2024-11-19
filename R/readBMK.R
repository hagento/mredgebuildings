#' Read heating systems sales in Austria from BMK
#'
#' @returns MagPIE object with yearly sales of heating systems in Germany
#'
#' @author Robin Hasse
#'
#' @importFrom utils read.csv
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate select matches
#' @importFrom magclass as.magpie collapseDim
#' @export

readBMK <- function() {
  read.csv("InnovativeEnergietechnologienInOesterreich2022_Abb29.txt",
           sep = "\t", encoding = "UTF-8") %>%
    rename(technology = "X") %>%
    pivot_longer(matches("\\d{4}"), names_to = "period") %>%
    mutate(period = as.numeric(sub("^X", "", .data[["period"]])),
           technology = sub(" +$", "", .data[["technology"]])) %>%
    as.magpie(temporal = "period", datacol = "value") %>%
    collapseDim(1)
}
