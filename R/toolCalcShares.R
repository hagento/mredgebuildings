#' Calculate Shares
#'
#' @param df data.frame
#' @param colShare column to be grouped over
#'
#' @returns data.frame with calculated shares for colShare
#'
#' @author Hagen Tockhorn
#'
#' @importFrom dplyr group_by across all_of mutate ungroup

toolCalcShares <- function(df, colShare) {
  df <- df %>%
    group_by(across(-all_of(c(colShare, "value")))) %>%
    mutate(value = replace_na(.data[["value"]], 0),
           value = proportions(.data[["value"]])) %>%
    ungroup()

  return(df)
}
