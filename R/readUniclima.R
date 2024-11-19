#' Read heating systems sales in France from Uniclima
#'
#' @returns magpie object
#'
#' @author Robin Hasse
#'
#' @importFrom utils read.csv
#' @importFrom dplyr %>% matches mutate select
#' @importFrom tidyr pivot_longer
#' @importFrom magclass as.magpie collapseDim
#' @export

readUniclima <- function() {
  read.csv("data.csv") %>%
    pivot_longer(matches("\\d{4}"), names_to = "period") %>%
    mutate(period = as.numeric(sub("^X", "", .data[["period"]]))) %>%
    select(-"comment") %>%
    as.magpie() %>%
    collapseDim(1)
}
