#' read EU28 building stock from Hotmaps
#'
#' Consistent building stock of EU28 residential and commercial sector. The data
#' presents a snapshot, I guess in 2016 but this is not clear to me yet.
#'
#' There are a few duplicate rows. Due to inconsistent values we average across
#' duplicates.
#'
#' @author Robin Hasse
#'
#' @importFrom utils read.csv
#' @importFrom dplyr %>% .data mutate select group_by summarise
#' @importFrom tidyr unite
#' @importFrom quitte as.quitte
#' @importFrom magclass as.magpie
#' @export

readHotmaps <- function() {

  data <- read.csv("building_stock.csv", sep = "|", encoding = "UTF-8") %>%
    unite("variable", "topic", "feature", "type", "detail", sep = "|") %>%
    unite("building", "sector", "subsector", "btype", sep = "|") %>%
    select(region = "country_code",
           "bage",
           "building",
           "variable",
           "value",
           "unit") %>%
    mutate(period = 2016,
           variable = sub(" *\\[.*\\]", "", .data[["variable"]]),
           variable = iconv(.data[["variable"]], "UTF-8", "ASCII", " "),
           variable = sub("\\|+$", "", .data[["variable"]]),
           variable = gsub("\\|{2,}", "|", .data[["variable"]]),
           variable = gsub(" {2,}", " ", .data[["variable"]]),
           variable = trimws(.data[["variable"]]),
           unit = sub("<c2><b2>", "2",
                      iconv(.data[["unit"]], "UTF-8", "ASCII", "byte"))) %>%
    # duplicate rows are averaged
    group_by(across(-all_of("value"))) %>%
    summarise(value = mean(.data[["value"]]), .groups = "drop") %>%
    as.quitte() %>%
    as.magpie()

  return(data)
}
