#' Temporary reading of HDD and CDD data
#'
#' Heating and cooling degree days calculated with \code{mredgebuildings::calcHDDCDD}.
#' Calculating this data currently requires manual help. We therefore read the
#' result as a source here. This becomes obsolete once the primary calculation
#' is stable.
#'
#' @author Robin Hasse
#'
#' @importFrom magclass as.magpie
#' @importFrom dplyr %>% .data group_by summarise ungroup
#' @export

readHDDCDDtemp <- function() {
  c("hddcddBAITnoCC_std2K.csv",
    "hddcddBAITssp2_std2K.csv") %>%
    lapply(read.csv) %>%
    do.call(what = rbind) %>%
    group_by(across(-all_of("value"))) %>%
    summarise(value = mean(.data[["value"]]), .groups = "drop") %>%
    ungroup() %>%
    as.magpie(temporal = "period", spatial = "region", datacol = "value")
}
