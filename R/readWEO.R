#' Read WEO enduse shares w.r.t. to global final energy demand
#'
#' @param subtype variable to define data subset ("Buildings", "Transport", "Industry")
#'
#' @author Hagen Tockhorn
#'
#' @importFrom tidyr gather
#' @importFrom dplyr mutate select filter
#' @importFrom magclass as.magpie
#' @importFrom quitte as.quitte
#' @importFrom readxl read_xlsx

readWEO <- function(subtype) {

  file <- "WEO_10.2.xlsx"

  # Read Data
  data <- read_xlsx(file, sheet = 2)

  # Process Data
  data <- data %>%
    filter(.data[["sector"]] == subtype) %>%
    gather("region", "value", 3:12) %>%
    mutate(period = 2014,
           region = gsub("\\.", " ", .data[["region"]]),
           enduse = tolower(.data[["enduse"]]),
           enduse = gsub(" ", "_", .data[["enduse"]])) %>%
    select("region", "period", "enduse", "value")

  data <- data %>%
    as.quitte() %>%
    as.magpie()

  return(data)
}
