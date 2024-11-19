#' Read heating systems stock and sales in Poland from SPIUG
#'
#' @returns magpie object
#'
#' @author Robin Hasse
#'
#' @importFrom utils read.csv
#' @importFrom dplyr %>% matches mutate
#' @importFrom tidyr pivot_longer
#' @importFrom magclass as.magpie collapseDim
#' @export

readSPIUG <- function() {

  data <- read.csv("heatingDeviceSales.csv") %>%
    pivot_longer(matches("\\d{4}"), names_to = "period") %>%
    mutate(period = as.numeric(sub("^X", "", .data[["period"]])))

  colnames(data)[1] <- "variable"

  as.magpie(data, temporal = "period") %>%
    collapseDim(1)
}
