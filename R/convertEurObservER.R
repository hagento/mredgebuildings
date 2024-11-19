#' convertEurObservER
#'
#' Each report of the Heat Pumps Barometer covers the two years before. There
#' numerous cases, where reports contradict in overlapping values or there
#' are huge jumps (probably due to changes in the accounting). This function
#' creates smooth trajectories back-casting from the most recent data. If data
#' overlaps, earlier data is scaled to match more recent data. In case of gaps,
#' the trajectories are back-casted maintaining relative growth within each
#' report. Growth rates between reports are interpolated. The Method yields
#' mostly smooth trajectories always hitting the mos recent reported numbers.
#'
#' @param x raw data
#' @param subtype Eurostat code of data set
#' @returns MAgPie object with converted data
#'
#' @author Robin Hasse
#'
#' @importFrom magclass getItems getItems<- getSets<- getSets
#' @importFrom madrat toolCountry2isocode
#' @export

convertEurObservER <- function(x, subtype) {


  backcast <- function(ts, key, growth) {

    # Check that there are only unique values in disregarded columns
    nUnique <- lapply(
      apply(select(ts, -any_of(c("period", "reportperiod", "value"))), 2, unique),
      length
    )
    if (any(nUnique > 1)) {
      stop(paste(names(nUnique)[which(nUnique > 1)], collapse = ", "),
           "have more than one unique value. Group by more dimensions.")
    }

    growthTs <- semi_join(growth, ts, by = c("region", "variable"))
    periods <- max(ts$period):min(ts$period)
    values <- c()
    v <- NA
    tRepsPrev <- c()

    for (i in seq_along(periods)) {
      t <- periods[i]
      tsStep <- filter(ts, .data[["period"]] == t)
      growthStep <- getElement(filter(growthTs, .data[["period"]] == t), "value")
      tReps <- as.numeric(as.character(getElement(tsStep, "reportperiod")))
      if (length(tReps) == 0) tReps <- NA

      if (i == 1 || max(tReps) %in% tRepsPrev) {
        v <- tsStep %>%
          filter(.data[["reportperiod"]] == max(tReps)) %>%
          getElement("value")
      } else {
        v <- v / (growthStep + 1)
      }

      values <- append(values, v)

      # scale overlapping reports to match
      ts <- ts %>%
        group_by(.data[["reportperiod"]]) %>%
        mutate(value = .data[["value"]] *
                 ifelse(.data[["reportperiod"]] %in% tReps,
                        (v / .data[["value"]][.data[["period"]] == t]),
                        1)) %>%
        ungroup()
      tRepsPrev <- tReps
    }

    # return backcasted time series
    ts %>%
      select(-"reportperiod", -"period", -"value", -"region", -"variable") %>%
      unique() %>%
      full_join(data.frame(period = periods, value = values),
                by = character())
  }


  data <- as.quitte(x) %>%
    filter(!is.na(.data[["value"]]))

  growth <- data %>%
    group_by(across(all_of(c("region", "reportperiod", "variable")))) %>%
    arrange(.data[["period"]]) %>%
    summarise(value = diff(.data[["value"]]) / diff(.data[["period"]]) / head(.data[["value"]], -1),
              period = head(.data[["period"]], -1),
              .groups = "drop") %>%
    select(-"reportperiod") %>%
    interpolate_missing_periods(min(data$period):max(data$period))

  # backcast
  data <- data %>%
    filter(.data[["variable"]] != "Total") %>%
    group_by(across(all_of(c("region", "variable")))) %>%
    group_modify(backcast, growth = growth, .keep = TRUE)

  # recalculate total
  data <- data %>%
    filter(.data[["variable"]] %in% c("ASHP", "GSHP")) %>%
    group_by(across(all_of(c("region", "period")))) %>%
    summarise(value = sum(.data[["value"]]),
              variable = "Total",
              .groups = "drop") %>%
    rbind(select(data, -"model", -"scenario", -"unit"))

  # convert to magpie object
  data <- data %>%
    as.quitte() %>%
    as.magpie()

  # rename regions: ISO2 -> ISO3
  getItems(data, 1) <- toolCountry2isocode(getItems(data, 1))

  # fill missing regions with NA
  data <- toolCountryFill(data, verbosity = 2)

  return(data)
}
