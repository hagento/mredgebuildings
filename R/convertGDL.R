#' convertGDL
#'
#' Drop subnational data but keep other subgroups (i.e. rural/urban, income
#' quartiles)
#'
#' @param subtype <database>.<variable>
#' @param x MAgPIE object with data from EU Buildings Database
#' @return clean MAgPIE object with data from EU Buildings Database
#'
#' @author Robin Hasse
#'
#' @importFrom magclass getItems getItems<- add_dimension getSets<-
#' @importFrom madrat toolCountryFill
#' @importFrom dplyr %>% filter mutate .data
#' @importFrom quitte as.quitte replace_column
#'
#' @export

convertGDL <- function(x, subtype) {

  data <- x

  # remove subnational regions
  subnationalRegs <- grep("^[A-Z]{3}r[0-9]*$", getItems(data, 1), value = TRUE)
  data <- data[subnationalRegs, , , invert = TRUE]

  # separate region (country) from subgroup
  data <- data %>%
    as.quitte() %>%
    filter(!is.na(.data[["value"]])) %>%
    mutate(subgroup = gsub("^[A-Z]{3}", "", .data[["region"]]),
           subgroup = gsub("^u$", "tu", .data[["subgroup"]]),
           region = substr(.data[["region"]], 1, 3)) %>%
    as.quitte() %>%
    as.magpie()

  # fill missing regions with NA
  getItems(data, 1) <- gsub("KSV", "KOS", getItems(data, 1))
  data <- toolCountryFill(data, verbosity = 2)

  return(data)
}
