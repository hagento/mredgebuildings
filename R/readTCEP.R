#' Read TCEP data base
#'
#' @param subtype character, type of data
#' @returns magpie object
#'
#' @author Hagen Tockhorn, Robin Hasse
#'
#' @importFrom dplyr %>% all_of mutate
#' @importFrom tidyr pivot_longer
#' @importFrom magclass as.magpie
#' @importFrom readxl read_xlsx
#' @importFrom quitte as.quitte
#' @export


readTCEP <- function(subtype = "enduse") {

  file <- switch(subtype,
    floorspace = "TCEP2014_figure_01_41_modified.xlsx",
    enduse = "TCEP2014_figure_01_39_data.xlsx",
    stop("unkown subtype: ", subtype)
  )

  data <- read_xlsx(file) %>%
    pivot_longer(-all_of(c("region", "period")), names_to = "variable") %>%
    mutate(unit = switch(subtype, enduse = "EJ", floorspace = "Billion m2")) %>%
    as.quitte() %>%
    as.magpie()

  return(data)
}
