#' Read Weibull lifetime distribution parameters for buildings
#'
#' @source https://doi.org/10.1016/j.jclepro.2019.118658
#'
#' @param subtype character, building subsector (either 'residential' or
#'   'commercial')
#' @returns MagPIE object with Weibull lifetime parameters for different world
#'   regions
#'
#' @author Robin Hasse
#'
#' @importFrom utils read.csv
#' @importFrom dplyr %>%
#' @importFrom tidyr pivot_longer
#' @importFrom magclass as.magpie
#' @export

readDeetman2020 <- function(subtype = "residential") {

  switch(subtype,
    # see Table 3 in Deetman et al. 2020
    residential = read.csv("table3.csv") %>%
      pivot_longer(-"region", names_to = "variable"),
    # see Deetman et al. 2020 p. 5
    commercial = data.frame(region = "GLO",
                            variable = c("shape", "scale", "mean_lt"),
                            value = c(1.44, 49.6, 45)),
    stop("invalid subtype: ", subtype)
  ) %>%
    as.magpie(spatial = "region", datacol = "value")
}
