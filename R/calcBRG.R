#' Calculate total boiler sales from BRG
#'
#' BRG has extremely comprehensive (and extremely expensive) sales data. But one
#' can recover total sales of boilers for space heating from the EHPA report.
#'
#' @return MagPIE object with total boiler sales in European countries
#'
#' @author Robin Hasse
#'
#' @importFrom madrat readSource
#' @importFrom magclass dimSums mselect
#' @importFrom dplyr %>% .data filter
#' @export

calcBRG <- function() {

  hpSales <- readSource("EHPA", subtype = "hpSales")
  hpMarketShare <- readSource("EHPA", subtype = "hpMarketShare")

  heatingHP <- toolGetMapping("technologyMapping_EHPA.csv",
                              type = "sectoral", where = "mredgebuildings") %>%
    filter(.data[["heatingHP"]]) %>%
    getElement("hpType")

  hpSales <- hpSales %>%
    mselect(corrected = TRUE, hpType = heatingHP) %>%
    dimSums(na.rm = TRUE)

  totalBoilerSales <- round(hpSales / hpMarketShare)

  return(list(x = totalBoilerSales,
              weight = NULL,
              min = 0,
              description = "total boiler sales for space heating",
              unit = "units/year"))
}
