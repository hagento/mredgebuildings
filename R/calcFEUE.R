#' Historical Final and Useful Energy Demand
#'
#' Historic Final Energy (FE) Demand and FE-EU-Efficiencies are taken to
#' compute Useful Energy (UE) Demand.
#'
#' @author Hagen Tockhorn
#'
#' @export


calcFEUE <- function() {
  # FUNCTIONS ------------------------------------------------------------------

  # Aggregate to Enduse-Level
  sumDF <- function(df, variables, newname) {
    enduseSum <- df %>%
      filter(.data[["enduse"]] %in% variables) %>%
      group_by(across(-any_of(c("value", "enduse")))) %>%
      summarise(value = sum(.data[["value"]], na.rm = TRUE), .groups = "drop") %>%
      mutate(enduse = newname) %>%
      ungroup() %>%
      select(all_of(colnames(df)))
    df %>%
      filter(!(.data[["enduse"]] %in% variables)) %>%
      rbind(enduseSum)
  }


  # READ-IN DATA ---------------------------------------------------------------

  fe <- calcOutput("FEbyEUEC", aggregate = FALSE) %>%
    as.quitte()


  efficiencies <- calcOutput("FEUEefficiencies", aggregate = FALSE) %>%
    as.quitte()


  # PROCESS DATA ---------------------------------------------------------------

  efficiencies <- rename(efficiencies, efficiency = "value")

  ue <- fe %>%
    sumDF(c("appliances", "lightning"), "appliances_light") %>%
    spread("unit", "value") %>%
    left_join(efficiencies,
              by = c("model", "scenario", "region", "variable",
                     "period", "carrier", "enduse")) %>%
    mutate(ue = .data[["fe"]] * .data[["efficiency"]]) %>%
    gather("unit", "value", "fe", "ue") %>%
    select(-"efficiency")


  # OUTPUT ---------------------------------------------------------------------

  ue <- ue %>%
    select("region", "period", "carrier", "enduse", "unit", "value") %>%
    as.quitte() %>%
    as.magpie()

  return(list(
    x = ue,
    weight = NULL,
    unit = "EJ",
    min = 0,
    description = "Final and Useful Energy Demand scaled with FE-UE-efficiencies"
  ))


}
