#' Historic Final Energy Demand by Carrier and Enduse
#'
#' Takes the historic final energy demand by carriers from IEA and disaggregates
#' it into different end uses.
#'
#' @returns data.frame with historic energy demands
#'
#' @author Hagen Tockhorn
#'
#' @importFrom dplyr mutate inner_join right_join left_join as_tibble filter
#'   select %>% .data
#' @importFrom madrat toolCountryFill
#' @importFrom quitte as.quitte
#' @importFrom magclass as.magpie
#'
#' @export


getFEbyEUEC <- function() {
  #---Read-in FE Data
  ieaIO <- calcOutput("IOEdgeBuildings", subtype = "output_EDGE_buildings",
                      aggregate = FALSE) %>%
    as.quitte(na.rm = TRUE)


  #---Read in Shares
  shares <- calcOutput("Shares", aggregate = FALSE) %>%
    as.quitte()


  #---Reduce the data frames dimensions to the minimal set
  commonRegionsPeriods <- Reduce(inner_join,
                                 list(unique(ieaIO[, c("region", "period")]),
                                      unique(shares[, c("region", "period")]))) %>%
    as_tibble()

  ieaIO <- commonRegionsPeriods %>%
    left_join(ieaIO, by = c("region", "period"))
  shares <- commonRegionsPeriods %>%
    left_join(shares, by = c("region", "period"))

  #---Calculate FE with shares
  ieaIO <- ieaIO %>%
    mutate(carrier = .data[["variable"]]) %>%
    select(-"model", -"scenario", -"variable", -"unit") %>%
    right_join(shares %>%
                 mutate(share = .data[["value"]]) %>%
                 select(-"model", -"scenario", -"variable", -"value"),
               by = c("region", "period", "carrier")) %>%
    mutate(value = .data[["share"]] * .data[["value"]],
           unit = "fe") %>%
    select(-"share")

  #---Pack Data
  data <- ieaIO %>%
    as.quitte() %>%
    as.magpie() %>%
    toolCountryFill(1, verbosity = 2)

  return(data)

}
