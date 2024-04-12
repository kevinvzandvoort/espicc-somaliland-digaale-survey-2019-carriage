#' set to adjust to deal with sparse data in pilot shipment
options(survey.lonely.psu = "adjust")

#############################################################################################################################
#' Section 1. Calculate carriage prevalence for these weights
#############################################################################################################################
table_data = rbind(
  data.table(variable_name = "Pneumococcal carriage",
    n = participant_data_design_svy$variables[!is.na(pneu_carr_final) & pneu_carr_final == 1, .N],
    N = participant_data_design_svy$variables[!is.na(pneu_carr_final), .N],
    type = "percentage") %>%
    .[, sample_mean := round(n/N * 100, 1)] %>%
    .[, c("pop_mean", "pop_confint_low", "pop_confint_high") :=
        svyMean2(~pneu_carr_final, participant_data_design_ps, multiply = 100, digits=1, na.rm=T) %>%
        .[, c("mean", "ci95_low", "ci95_high")]] %>%
    .[, c("popu5_mean", "popu5_confint_low", "popu5_confint_high") :=
        svyMean2(~pneu_carr_final, participant_data_u5_design_ps, multiply = 100, digits=1, na.rm=T) %>%
        .[, c("mean", "ci95_low", "ci95_high")]])
table_data[, heading := POSTSTRATIFICATION_STRATA]  

#' These are the new adjusted estimates
prevalence_data = svyMean2(~pneu_carr_final, participant_data_design_ps, na.rm=T, by="agegrp")
prevalence_data_bytype = svyMean2(~pneumosil_type, participant_data_design_ps, na.rm=T, by="agegrp") %>%
  .[, c("agegrp", "option", "mean")] %>%
  .[option != "S"] %>%
  .[, variable := option]

prevalence_data[, stratum := POSTSTRATIFICATION_STRATA]
prevalence_data_bytype[, stratum := POSTSTRATIFICATION_STRATA]

postratification_sensitivity_data[[POSTSTRATIFICATION_STRATA]] = list(table_data = table_data,
                                                                      plot_data = list(
                                                                        prevalence_data = prevalence_data,
                                                                        prevalence_data_bytype = prevalence_data_bytype))

#' return to default
options(survey.lonely.psu = "fail")