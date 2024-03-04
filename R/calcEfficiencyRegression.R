#' Calculation Regression Parameters for FE-UE Efficiency Projection
#'
#' Calculate the regression parameters for the projection of final (FE)
#' to useful (UE) energy conversion projection for all combinations of
#' energy enduses and carriers. The regression parameters correspond to an
#' asymptotic regression model and encompass the parameters Asym, R0 and lrc.
#' The are determined using a nonlinear least-squares regression.
#'
#' This approach closely follows the model by De Stercke et al. which is
#' mainly driven by GDP per Capita.
#'
#' @references De Stercke, S. (2014). Dynamics of Energy Systems: A Useful
#' Perspective (Interim Report, p. 68). IIASA.
#' http://pure.iiasa.ac.at/id/eprint/11254/1/IR-14-013.pdf
#'
#' @author Hagen Tockhorn
#'
#' @importFrom quitte aggregate_map
#' @importFrom stats as.formula
#'
#' @export


calcEfficiencyRegression <- function() {
  # FUNCTIONS ------------------------------------------------------------------

  # Extrapolate historic FE-UE Efficiencies from Fit Function
  calcPars <- function(df, var) {
    # Prepare Historic Data
    dataHist <- df %>%
      removeColNa() %>%
      filter(.data[["variable"]] == var) %>%
      spread("variable", "value") %>%
      na.omit()

    # Replace 0's with small Values to avoid Inf Issues
    dataHist[dataHist[var] == 0, var] <-
      min(dataHist[dataHist[var] != 0, var]) / 10

    # Create Estimation Object for Non-Linear Model
    estimate <- nls(as.formula(paste(var, "~ SSasymp(gdppop, Asym, R0, lrc)")),
                    dataHist)

    return(estimate$m$getPars())
  }



  # READ-IN DATA ---------------------------------------------------------------

  pfu <- calcOutput("PFUDB", aggregate = FALSE) %>%
    as.quitte()

  gdppop <- calcOutput("GDPPop", aggregate = FALSE) %>%
    as.quitte()

  # Get Mapping (ISO<->PFU)
  regionmapping <- toolGetMapping("pfu_regionmapping.csv", type = "regional")

  # Get Population Data
  pop <- calcOutput("PopulationPast", aggregate = FALSE) %>%
    as.quitte()



  # PARAMETERS -----------------------------------------------------------------

  # Minimum Requirement to be considered
  minEfficiency <- 0.05



  # PROCESS DATA ---------------------------------------------------------------

  #--- Calculate existing FE-EU efficiencies

  # Aggregate PFU Data to PFU Country Code
  data <- pfu %>%
    mutate(value = ifelse(is.na(.data[["value"]]), 0, .data[["value"]])) %>%
    unite("variable", "enduse", "carrier", sep = ".") %>%
    quitte::aggregate_map(
      mapping = regionmapping[!is.na(regionmapping$PFUDB), c("iso", "PFUDB")],
      by = c("region" = "iso"),
      forceAggregation = TRUE)

  # Aggregate GDPpop to PFU Country Code
  gdppop <- gdppop %>%
    quitte::aggregate_map(
      mapping = regionmapping[!is.na(regionmapping$PFUDB), c("iso", "PFUDB")],
      by = c("region" = "iso"),
      forceAggregation = TRUE,
      weights = pop %>%
        select("region", "period", "value") %>%
        rename(weight = "value"),
      weight_item_col = "region",
      weight_val_col = "weight") %>%
    select(-"model", -"scenario", -"unit", -"variable")


  # Combine with GDP per Cap
  data <- data %>%
    left_join(gdppop %>%
                rename(gdppop = "value"),
              by = c("region", "period"))

  # Calculate Efficiencies as Regression-Input
  data <- data %>%
    select(-"model", -"scenario") %>%
    spread(.data[["unit"]], .data[["value"]]) %>%
    mutate(value = .data[["ue"]] / .data[["fe"]]) %>%
    select(-"fe", -"ue")

  # Filter out unrealistic Efficiencies
  data <- filter(data, .data[["value"]] > minEfficiency)

  vars <- data %>%
    group_by(.data[["variable"]]) %>%
    filter(.data[["variable"]] != "gdppop",
           !all(is.na(.data[["value"]]))) %>%
    getElement("variable") %>%
    unique()


  #--- Calculate Parameters
  parsFull <- data.frame()

  for (var in vars) {
    pars <- calcPars(data, var)
    parsFull <- as.data.frame(do.call(cbind, as.list(pars))) %>%
      mutate(variable = var) %>%
      rbind(parsFull)
  }



  # OUTPUT ---------------------------------------------------------------------

  # Trim Dataframe
  parsFull <- parsFull %>%
    separate("variable", c("enduse", "carrier"), sep = "\\.") %>%
    mutate(region = "GLO") %>%
    select("region", "carrier", "enduse", "Asym", "R0", "lrc")


  return(list(
    x = as.magpie(parsFull),
    weight = NULL,
    description = "Regression Parameter for FE-UE-Efficiency Projection",
    unit = ""
  ))

}
