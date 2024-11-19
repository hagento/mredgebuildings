#' Read EHI Heating Market Report
#'
#' @param subtype Character, report to read. if NULL, the most recent report is
#'   used for every data point.
#' @returns magpie object
#'
#' @author Robin Hasse
#'
#' @importFrom utils read.csv
#' @importFrom dplyr %>% mutate filter select across all_of row_number arrange n
#' @importFrom magclass as.magpie
#' @export

readEHI <- function(subtype = NULL) {

  # names of files must be orderable
  files <- c(`2021` = "heatingMarketReport2021.csv",
             `2023` = "heatingMarketReport2023.csv")

  data <- if (is.null(subtype)) {
    # consider all files and always use the most recent information
    do.call(rbind, lapply(names(files), function(fileName) {
      read.csv(files[[fileName]]) %>%
        mutate(file = fileName)
    })) %>%
      group_by(across(-all_of(c("file", "value")))) %>%
      arrange(.data[["file"]]) %>%
      filter(row_number() == n()) %>%
      ungroup() %>%
      select(-"file") %>%
      as.magpie(spatial = "region", temporal = "period", datacol = "value")
  } else if (subtype %in% names(files)) {
    # read only given report
    read.csv(files[[subtype]]) %>%
      as.magpie(spatial = "region", temporal = "period", datacol = "value")
  } else if (grepl("^consumerStudy2021\\.", subtype)) {
    file <- paste0(sub("^consumerStudy2021\\.", "", subtype), ".csv")
    replacementMatrix <- grepl("^table2-[2-6]\\.csv", file)
    data <- read.csv(file.path("consumerStudy2021", file),
                     encoding = "UTF-8", check.names = FALSE)
    colnames(data) <- ifelse(colnames(data) == "", seq_along(data), colnames(data))
    if (replacementMatrix) {
      # remove total of old shares
      data[[2]] <- NULL
    }
    data <- data %>%
      pivot_longer(-1) %>%
      mutate(value = as.numeric(sub("\\%$", "", .data[["value"]])) / 100)
    if (replacementMatrix) {
      colnames(data)[1:2] <- c("old", "new")
      data <- as.magpie(data, datacol = "value") %>%
        collapseDim(c(1, 2))
    } else {
      colnames(data)[1:2] <- c("age", "region")
      data <- as.magpie(data, spatial = "region", datacol = "value") %>%
        collapseDim(2)
    }
  } else {
    stop("Invalid subtype: '", subtype, "'. Available subtypes: NULL, ",
         paste(names(files), collapse = ", "))
  }

  return(data)
}
