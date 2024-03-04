#' Calculation and Projection of Final to Useful Energy
#'
#' Calculate Efficiencies of Final (FE) to Useful (UE) Energy Conversion for all
#' combinations of Energy Carriers and Enduses.
#' The efficiency projections are based on a model by De Stercke et al. which is
#' mainly driven by GDP per Capita. It describes an S-shaped curve approaching
#' assumed efficiency levels. The parameters of that curve are derived by a
#' regression with observations of IEA data.
#'
#' @param gasBioEquality Determines if Natgas and Biomod share the same efficiencies
#'
#' @references De Stercke, S. (2014). Dynamics of Energy Systems: A Useful
#' Perspective (Interim Report, p. 68). IIASA.
#' http://pure.iiasa.ac.at/id/eprint/11254/1/IR-14-013.pdf
#'
#' @author Hagen Tockhorn
#'
#' @importFrom stats SSasymp
#'
#' @export


calcFEUEefficiencies <- function(gasBioEquality = TRUE) {
  # READ-IN DATA ---------------------------------------------------------------

  pfu <- calcOutput("PFUDB", aggregate = FALSE) %>%
    as.quitte()

  gdppop <- calcOutput("GDPPop", aggregate = FALSE) %>%
    as.quitte()

  # parsCorr <- read.csv("correct_efficiencies.csv", sep = ";")
  parsCorr <- toolGetMapping("correct_efficiencies.csv", type = "sectoral")


  # PARAMETERS -----------------------------------------------------------------

  eqEffs <- c("water_heating.natgas" = "water_heating.biomod",
              "space_heating.natgas" = "space_heating.biomod",
              "cooking.natgas" = "cooking.biomod")


  # PROCESS DATA ---------------------------------------------------------------

  # Combine with GDP per Cap for Period 1990-2020
  data <- pfu %>%
    interpolate_missing_periods(period = seq(1990, 2020)) %>%
    mutate(value = ifelse(is.na(.data[["value"]]), 0, .data[["value"]])) %>%
    left_join(gdppop %>%
                select(-"model", -"scenario", -"unit", -"variable") %>%
                rename(gdppop = "value"),
              by = c("region", "period")) %>%
    select(-"model", -"scenario", -"variable")


  #--- Calculate Efficiency Estimates

  # Regression Parameters for enduse.carrier Combinations
  regPars <- calcOutput("EfficiencyRegression", aggregate = FALSE) %>%
    as.quitte() %>%
    spread(key = "variable", value = "value") %>%
    select(-"model", -"scenario", -"region", -"unit", -"period")


  # Tediously transform the corrected parameters ... must be a better solution
  parsCorr <- parsCorr %>%
    mutate("space_cooling.elec" = gsub(" ", "", .data[["space_cooling.elec"]]),
           "space_cooling.elec" = gsub("\\,", "\\.", .data[["space_cooling.elec"]]),
           "space_cooling.elec" = as.numeric(.data[["space_cooling.elec"]])) %>%
    gather(key = "variable", value = "value", "space_cooling.elec", "water_heating.elec", "space_heating.elec") %>%
    spread(key = "parameter", value = "value") %>%
    separate(col = "variable", into = c("enduse", "carrier"), sep = "\\.") %>%
    select(-"phi3") %>%
    rename(AsymCorr = "Asym", lrcCorr = "lrc", R0Corr = "R0")



  # Correct Regression Parameters for electrical Heat Transfer Technologies
  regPars <- regPars %>%
    left_join(parsCorr, by = c("carrier", "enduse")) %>%
    mutate(Asym = ifelse(is.na(.data[["AsymCorr"]]), .data[["Asym"]], .data[["AsymCorr"]]),
           R0 = ifelse(is.na(.data[["R0Corr"]]), .data[["R0"]], .data[["R0Corr"]]),
           lrc = ifelse(is.na(.data[["lrcCorr"]]), .data[["lrc"]], .data[["lrcCorr"]])) %>%
    select(-"AsymCorr", -"R0Corr", -"lrcCorr")


  # Assign Parameters to carrier-enduse Combination
  dataHist <- data %>%
    removeColNa() %>%
    left_join(regPars, by = c("carrier", "enduse")) %>%
    na.omit()


  # Predict Historic Efficiencies with Non-Linear Model
  dataHist <- dataHist %>%
    spread(key = "unit", value = "value") %>%
    mutate(efficiency = .data[["ue"]] / .data[["fe"]]) %>%
    group_by(across(all_of(c("carrier", "enduse")))) %>%
    mutate(pred = SSasymp(.data[["gdppop"]], .data[["Asym"]], .data[["R0"]], .data[["lrc"]])) %>%
    ungroup() %>%
    select(-"Asym", -"R0", -"lrc", -"fe", -"ue")


  #--- Match Region-Specific Curves to fill non-existing Data Points

  # Create Correction Factor to adjust Projections
  corrFactors <- dataHist %>%
    mutate(factor = .data[["efficiency"]] / .data[["pred"]]) %>%
    select(-"gdppop", -"efficiency", -"pred") %>%
    interpolate_missing_periods(value = "factor", expand.values = TRUE)


  # NOTE: Countries missing the entire period range for a EC-EU-combination
  #       will be filled-up w/ non-corrected efficiency projections.
  dataHist <- dataHist %>%
    left_join(corrFactors, by = c("region", "period", "enduse", "carrier")) %>%
    mutate(value = ifelse(is.na(.data[["factor"]]),
                          .data[["pred"]],
                          ifelse(is.infinite(.data[["factor"]]),
                                 .data[["pred"]],
                                 ifelse(is.na(.data[["efficiency"]]),
                                        .data[["pred"]] * .data[["factor"]],
                                        .data[["efficiency"]])
                                 )
                          )
           )


  # Trim Dataframe
  efficiencies <- dataHist %>%
    select(-"gdppop", -"efficiency", -"pred", -"factor")


  #--- Corrections

  # Biomod Efficiency identical to Natgas
  if (gasBioEquality) {
    gasEffs <- efficiencies %>%
      unite("variable", "enduse", "carrier", sep = ".") %>%
      filter(.data[["variable"]] %in% names(eqEffs))

    for (gasVar in names(eqEffs)) {
      bioEffs <- gasEffs %>%
        filter(.data[["variable"]] == gasVar) %>%
        mutate(variable = eqEffs[gasVar][[1]]) %>%
        separate("variable", into = c("enduse", "carrier"), sep = "\\.")

      efficiencies <- rbind(efficiencies, bioEffs)
    }
  }


  # OUTPUT ---------------------------------------------------------------------

  # Weights = FE?

  efficiencies <- efficiencies %>%
    as.quitte() %>%
    as.magpie()

  return(list(
    x = efficiencies,
    weights = NULL,
    unit = "",
    min = 0,
    description = "Historical Conversion Efficiencies from FE to UE"
  ))






}
