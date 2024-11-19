#' Read heat pump sales data from EHPA
#'
#' This data was provided upon request. It is free for non-profit organisations
#' but we are not allowed to share it.
#'
#' @param subtype character, type of data
#' @returns MagPIE object with HP sales data for EU member states
#'
#' @author Robin Hasse
#'
#' @importFrom utils read.csv
#' @importFrom dplyr rename mutate matches
#' @importFrom tidyr pivot_longer
#' @importFrom magclass as.magpie collapseDim
#' @export

readEHPA <- function(subtype = "hpSales") {

  switch(subtype,
    hpSales = {

      file <- "EHPA_sales_full sheet from stats tool_Dec 2023.csv"

      # the file has weired quotations
      text <- readLines(file)
      text <- gsub("^\"|\"$", "", text)
      text <- gsub("\"\"", "\"", text)

      data <- read.csv(text = text) %>%
        rename(hpType = "HP.type") %>%
        pivot_longer(matches("^sales"), names_to = "corrected") %>%
        mutate(corrected = !grepl("no.correction", .data[["corrected"]])) %>%
        as.magpie(spatial = "country", temporal = "year", datacol = "value")
    },
    hpMarketShare = {

      file <- "EHPA_market_report_2023_Final full report_Table3.5-1.txt"
      data <- read.csv(file, sep = " ", comment.char = "#") %>%
        rename(region = "X") %>%
        pivot_longer(matches("X\\d{4}"), names_to = "period") %>%
        mutate(period = as.numeric(sub("^X", "", .data[["period"]])),
               value = as.numeric(sub("%$", "", .data[["value"]])) / 100) %>%
        as.magpie(spatial = "region", temporal = "period", datacol = "value") %>%
        collapseDim()
    },
    stop("Unkown subtype")
  )

  return(data)
}
