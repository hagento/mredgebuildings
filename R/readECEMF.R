#' Read ECEMF data
#'
#' Project data: European Climate and Energy Modelling Forum
#'
#' - FEPrices: price components w/o VAT in EUR/MWh
#'
#' @author Robin Hasse
#'
#' @param subtype character, type of data
#'
#' @importFrom openxlsx read.xlsx
#' @importFrom tidyr pivot_longer replace_na
#' @importFrom dplyr select mutate matches filter
#' @export

readECEMF <- function(subtype = "FEPrices") {

  # raw data
  ecemf <- read.xlsx("FE Prices_2022-04-14_v2.xlsx", "Price markups",
                     fillMergedCells = TRUE)
  switch(subtype,
    FEPrices = {
      # only read price components
      data <- ecemf %>%
        select(carrier = "X1", sector = "X3", component = "X4", unit = "X7",
               matches("\\d{4}")) %>%
        mutate(carrier = ifelse(grepl("average of model ensemble|from CO2 price path",
                                      .data[["carrier"]]),
                                NA,
                                .data[["carrier"]]),
               pos = cumsum(!is.na(.data[["carrier"]]))) %>%
        filter(.data[["pos"]] > 0) %>%
        mutate(carrier = unique(.data[["carrier"]][!is.na(.data[["carrier"]])])[.data[["pos"]]]) %>%
        filter(.data[["unit"]] == "EUR/MWh") %>%
        select(-"pos", -"unit") %>%
        mutate(component = ifelse(grepl("^Price (PE|Pe|SE)", .data[["component"]]),
                                  "energy",
                                  .data[["component"]])) %>%
        pivot_longer(matches("\\d{4}"), names_to = "period", values_to = "value") %>%
        mutate(period = as.numeric(.data[["period"]]),
               sector = replace_na(.data[["sector"]], "Supply"))
    },
    emissionFactor = {
      data <- ecemf %>%
        select(carrier = "X1", variable = "X4", value = "X5") %>%
        filter(.data[["variable"]] == "Emission factor") %>%
        mutate(value = as.numeric(.data[["value"]]))
    },
    stop("Invalid subtype: ", subtype)
  )

  # convert to MAgPIE object
  data <- data %>%
    as.magpie() %>%
    collapseDim()

  return(data)

}
