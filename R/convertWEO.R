#' Convert WEO enduse shares w.r.t. to global final energy demand
#'
#' @param x readWEO object
#'
#' @author Hagen Tockhorn
#'
#' @importFrom madrat toolGetMapping
#' @importFrom quitte aggregate_map

convertWEO <- function(x) {

  # PARAMETER-------------------------------------------------------------------

  # small correction factor for 0 entries
  eps <- 1e-5

  # READ DATA-------------------------------------------------------------------

  data <- as.quitte(x)

  # region mapping
  regmapping <- toolGetMapping(name  = "regionmappingWEO.csv",
                               type  = "regional",
                               where = "mredgebuildings")

  # fe weights
  fe <- calcOutput("ShareETP",
                   subtype = "enduse",
                   aggregate = FALSE,
                   feOnly = TRUE) %>%
    as.quitte()

  # correction values
  correction <- toolGetMapping(name  = "correct_WEO.csv",
                               type  = "sectoral",
                               where = "mredgebuildings")



  # CORRECTIONS-----------------------------------------------------------------

  data <- data %>%
    select("region", "period", "enduse", "value") %>%
    rename(variable = "enduse") %>%
    mutate(unit = NA)

  dataCorr <- data %>%
    left_join(correction, by = c("region", "variable")) %>%
    mutate(value = ifelse(is.na(.data[["correction"]]),
                          .data[["value"]],
                          .data[["value"]] * .data[["correction"]])) %>%
    select(-"correction")


  # PROCESS DATA----------------------------------------------------------------

  fe <- fe %>%
    select("region", "period", "enduse", "value") %>%
    rename("variable" = "enduse")


  dataDisagg <- dataCorr %>%
    aggregate_map(mapping = regmapping,
                  by = c(region = "RegionCode"),
                  weights = fe %>%
                    rename(weight_val_col = "value"),
                  weight_item_col = "region",
                  forceAggregation = TRUE) %>%
    rename(enduse = "variable")


  data <- dataDisagg %>%
    select("region", "period", "enduse", "value")

  # calculate regional shares
  data <- data %>%
    group_by(across(all_of(c("region", "period")))) %>%
    mutate(value = proportions(.data[["value"]])) %>%
    ungroup()


  # OUTPUT----------------------------------------------------------------------

  data <- data %>%
    as.quitte() %>%
    as.magpie %>%
    toolCountryFill()

  return(data)
}
