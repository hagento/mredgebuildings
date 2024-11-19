#' Introduce refrigerators as appliances
#'
#' @param df data.frame
#' @param mapping original EDGE-B region mapping
#' @param fridgeShare regional shares of refrigerators
#' @param feOnly are absolute values (TRUE) or shares (FALSE) considered?
#' @param shareOf column name to calculate share of - must be given if feOnly = FALSE
#'
#' @author Hagen Tockhorn
#'
#' @importFrom dplyr filter left_join select

toolAddThermal <- function(df, mapping, fridgeShare, feOnly = TRUE, shareOf = NULL) {
  df <- df %>%
    filter(.data[["enduse"]] != "lighting") %>%
    left_join(mapping %>%
                select(-"RegionCodeEUR", -"RegionCodeEUR_ETP", -"X") %>%
                rename(region = "CountryCode") %>%
                left_join(fridgeShare, by = "RegionCode") %>%
                select(-"RegionCode"),
              by = "region") %>%
    mutate(value = ifelse(.data[["enduse"]] != "appliances",
                          .data[["value"]],
                          .data[["value"]] * .data[["share"]]),
           enduse = ifelse(.data[["enduse"]] == "appliances",
                           "refrigerators",
                           as.character(.data[["enduse"]]))) %>%
    select(-"share")

  if (!feOnly) {
    if (is.character(shareOf) && (shareOf %in% colnames(df))) {
      df <- df %>%
        group_by(across(-all_of(c(shareOf, "value")))) %>%
        mutate(value = .data[["value"]] / sum(.data[["value"]], na.rm = TRUE)) %>%
        ungroup()
    } else {
      stop("Please give existing col name.")
    }

    return(df)
  }

  return(df)
}
