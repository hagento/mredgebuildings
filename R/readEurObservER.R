#' Read EurObservER
#'
#' PDF market reports on Renewables in Europe
#'
#' @source https://www.eurobserv-er.org/
#' @note I download the report and snip a png of the table. Then I extract the
#' data with [this online tool](https://extracttable.com/) and save it as a
#' tabular csv in the source folder. The Tool is not fully reliable. This
#' function throws an error if values don't sum up. You then have to correct the
#' csv file manually.
#'
#' @param subtype <report series>.<variable>
#' @returns magpie object
#'
#' @author Robin Hasse
#'
#' @importFrom utils read.csv capture.output
#' @importFrom dplyr %>% select mutate group_by filter rename .data
#' @importFrom tidyr gather
#' @importFrom quitte as.quitte
#' @importFrom magclass as.magpie
#' @export

readEurObservER <- function(subtype) {

  # PREPARE --------------------------------------------------------------------

  # split subtype
  report   <- gsub("\\..*$", "", subtype)

  # list report files
  if (!report %in% dir()) {
    stop("The subtype '", subtype, "' is invalid. Available reports: ",
         paste(dir(), collapse = ", "))
  }
  folder <- file.path(report)
  files <- list.files(folder, "\\.csv$", recursive = TRUE, full.names = TRUE)



  # PROCESS DATA ---------------------------------------------------------------

  switch(subtype,
    heatPumpsBarometer.stock = {

      # read and join all files
      files <- grep("(20(1[4-9]|2[0-9])_Table3|2013_Table[5-6])\\.csv$", files, value = TRUE)

      data <- do.call("rbind", lapply(files, function(file) {
        reportPeriod <- gsub(".*(\\d{4}).*", "\\1", file)
        wideTable <- read.csv(file, header = FALSE)

        if (reportPeriod == "2013") {
          years <- switch(gsub(".*Table|\\.csv", "", file), `5` = 2011, `6` = 2012)
          wideTable <- wideTable[, !grepl("[rR]enewable heat", wideTable[1, ])]
        } else {
          years <- grep("\\d{4}", wideTable[1, ], value = TRUE)
        }

        nDataCols <- (ncol(wideTable) - 1) / length(years)

        longTable <- do.call("rbind", lapply(seq_along(years), function(i) {
          subTable <- wideTable[ifelse(reportPeriod == "2013", 1, 2):nrow(wideTable),
                                c(1, ((i - 1) * nDataCols + 2):(i * nDataCols + 1))]
          colnames(subTable) <- subTable[1, ]
          subTable <- subTable[-1, ]
          subTable[["period"]] <- years[[i]]
          return(subTable)
        }))

        colnames(longTable)[1] <- "region"
        longTable <- longTable %>%
          gather("variable", "value", -any_of(c("region", "period"))) %>%
          mutate(reportPeriod = reportPeriod)

        return(longTable)
      }))

      # tidy data
      variableMatching <- c(
        `ASHP`                    = "ASHP",
        `Aerothermal HP`          = "ASHP",
        `Aerothermal heat pumps`  = "ASHP",
        `includ. air- water HP`   = "ASHP_airWater",
        `GSHP`                    = "GSHP",
        `Ground source heat pump` = "GSHP",
        `Geothermal HP`           = "GSHP",
        `Geothermal heat pumps`   = "GSHP",
        `Total PAC`               = "Total",
        `Total heat pumps`        = "Total",
        `Total HP in operation`   = "Total"
      )
      data <- data %>%
        mutate(value = ifelse(.data[["value"]] == "-", NA, .data[["value"]]),
               value = as.numeric(gsub("o|O", "0", gsub(" ", "", .data[["value"]]))),
               variable = variableMatching[.data[["variable"]]],
               period = as.numeric(.data[["period"]]),
               reportPeriod = as.numeric(.data[["reportPeriod"]]),
               region = iconv(.data[["region"]], "latin1", "ASCII", sub = ""),
               region = gsub("\\*", "", .data[["region"]]),
               region = gsub("Czechia", "Czech Republic", .data[["region"]]),
               region = gsub("^Netherland$", "Netherlands", .data[["region"]]),
               region = gsub("United(-K| k)ingdom", "United Kingdom", .data[["region"]]),
               region = gsub("Lituania", "Lithuania", .data[["region"]]),
               region = gsub("Hungarie", "Hungary", .data[["region"]]))
      sumIsWrong <- data %>%
        group_by(across(all_of(c("region", "period", "reportPeriod")))) %>%
        filter(abs(.data[["value"]][.data[["variable"]] == "Total"] -
                     .data[["value"]][.data[["variable"]] == "ASHP"] -
                     .data[["value"]][.data[["variable"]] == "GSHP"]) > 0) %>%
        select("region", "period", "reportPeriod") %>%
        unique()
      if (nrow(sumIsWrong) > 0) {
        stop("Total doesn't add up in the following cases. Check CSV files against PNGs.\n",
             paste(capture.output(print(sumIsWrong)), collapse = "\n"))
      }
    },
    heatPumpsBarometer.sales = {

      # read and join all files
      files <- c(grep("\\d{4}_Table(1|2)\\.csv$", files, value = TRUE))

      data <- do.call("rbind", lapply(files, function(file) {
        reportPeriod <- gsub(".*(\\d{4}).*", "\\1", file)
        wideTable <- read.csv(file, header = FALSE)

        if (reportPeriod == "2013") {
          years <- switch(gsub(".*Table|.csv", "", file), `5` = 2011, `6` = 2012)
          wideTable <- wideTable[, !grepl("[rR]enewable heat", wideTable[1, ])]
        } else {
          years <- grep("\\d{4}", wideTable[1, ], value = TRUE)
        }

        nDataCols <- (ncol(wideTable) - 1) / length(years)

        longTable <- do.call("rbind", lapply(seq_along(years), function(i) {
          subTable <- wideTable[ifelse(reportPeriod == "2013", 1, 2):nrow(wideTable),
                                c(1, ((i - 1) * nDataCols + 2):(i * nDataCols + 1))]
          colnames(subTable) <- subTable[1, ]
          subTable <- subTable[-1, ]
          subTable[["period"]] <- years[[i]]
          return(subTable)
        }))

        colnames(longTable)[1] <- "region"
        longTable <- longTable %>%
          gather("variable", "value", -any_of(c("region", "period"))) %>%
          mutate(reportPeriod = reportPeriod)

        return(longTable)
      }))

      # tidy data
      variableMatching <- c(
        `ASHP`                    = "ASHP",
        `Aerothermal HP`          = "ASHP",
        `Aerothermal heat pumps`  = "ASHP",
        `includ. air- water HP`   = "ASHP_airWater",
        `GSHP`                    = "GSHP",
        `Ground source heat pump` = "GSHP",
        `Geothermal HP`           = "GSHP",
        `Geothermal heat pumps`   = "GSHP",
        `Total PAC`               = "Total",
        `Total heat pumps`        = "Total",
        `Total HP in operation`   = "Total"
      )
      data <- data %>%
        mutate(value = ifelse(.data[["value"]] == "-", NA, .data[["value"]]),
               value = as.numeric(gsub("o|O", "0", gsub(" ", "", .data[["value"]]))),
               variable = variableMatching[.data[["variable"]]],
               period = as.numeric(.data[["period"]]),
               reportPeriod = as.numeric(.data[["reportPeriod"]]),
               region = iconv(.data[["region"]], "latin1", "ASCII", sub = ""),
               region = gsub("\\*", "", .data[["region"]]),
               region = gsub("Czechia", "Czech Republic", .data[["region"]]),
               region = gsub("^Netherland$", "Netherlands", .data[["region"]]),
               region = gsub("United(-K| k)ingdom", "United Kingdom", .data[["region"]]),
               region = gsub("Lituania", "Lithuania", .data[["region"]]),
               region = gsub("Hungarie", "Hungary", .data[["region"]]))
      sumIsWrong <- data %>%
        group_by(across(all_of(c("region", "period", "reportPeriod")))) %>%
        filter(abs(.data[["value"]][.data[["variable"]] == "Total"] -
                     .data[["value"]][.data[["variable"]] == "ASHP"] -
                     .data[["value"]][.data[["variable"]] == "GSHP"]) > 1) %>%
        select("region", "period", "reportPeriod") %>%
        unique()
      if (nrow(sumIsWrong) > 0) {
        stop("Total doesn't add up in the following cases. Check CSV files against PNGs.\n",
             paste(capture.output(print(sumIsWrong)), collapse = "\n"))
      }
    },
    stop("'", subtype, "' is an unsupported subtype.")
  )

  data <- data %>%
    as.quitte() %>%
    as.magpie()

  return(data)
}
