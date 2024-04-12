#' Create the different sections used in table 1

#############################################################################################################################
#' Section 1. Demographic characteristics of the sample
#############################################################################################################################
table1.1 = rbind(
  data.table(variable_name = "Median age (unweighted)",
    n = NA_integer_,
    N = participant_data_design_svy$variables[!is.na(participant_age_y), .N],
    sample_mean = svyquantile(~participant_age_y, participant_data_design_svy, c(0.5), na.rm = TRUE)[[1]][1,1],
    sample_confint_low = svyquantile(~participant_age_y, participant_data_design_svy, c(0.25), na.rm = TRUE)[[1]][1,1],
    sample_confint_high = svyquantile(~participant_age_y, participant_data_design_svy, c(0.75), na.rm = TRUE)[[1]][1,1],
    type = "IQR"),
  data.table(variable_name = "Percentage Female (unweighted)",
    n = participant_data_design_svy$variables[!is.na(participant_sex), sum(participant_sex == "female")],
    N = participant_data_design_svy$variables[!is.na(participant_sex), .N],
    sample_mean = svyMean2(~participant_sex, participant_data_design_svy, multiply = 100, digits=1, na.rm=T) %>%
      .[option == "female", mean],
    sample_confint_low = NA_real_,
    sample_confint_high = NA_real_,
    type = "percentage"),
  data.table(variable_name = "Median household size (unweighted)",
    n = NA_integer_,
    N = participant_data_design_svy$variables[!is.na(household_size), .N],
    sample_mean = svyquantile(~household_size, participant_data_design_svy, c(0.5), na.rm=T)[[1]][1,1],
    sample_confint_low = svyquantile(~household_size, participant_data_design_svy, c(0.25), na.rm=T)[[1]][1,1],
    sample_confint_high = svyquantile(~household_size, participant_data_design_svy, c(0.75), na.rm=T)[[1]][1,1],
    type = "IQR"),
  data.table(variable_name = "Median household members <5y (unweighted)",
    n = NA_integer_,
    N = participant_data_design_svy$variables[!is.na(household_size_u5), .N],
    sample_mean = svyquantile(~household_size_u5, participant_data_design_svy, c(0.5), na.rm=T)[[1]][1,1],
    sample_confint_low = svyquantile(~household_size_u5, participant_data_design_svy, c(0.25), na.rm=T)[[1]][1,1],
    sample_confint_high = svyquantile(~household_size_u5, participant_data_design_svy, c(0.75), na.rm=T)[[1]][1,1],
    type = "IQR"))
table1.1[, heading := "Sample characteristics"]

#############################################################################################################################
#' Section 2. Prevalence of potential risk factors
#' - Self-reported respiratory symptoms in the two weeks before swab collection
#' - contacts reported during the contact survey
#############################################################################################################################
table1.2 = rbind(  
  data.table(variable_name = "Respiratory symptoms",
    n = participant_data_design_ps$variables[!is.na(symptoms) & symptom_respiratory == "yes", .N],
    N = participant_data_design_ps$variables[!is.na(symptoms) & symptoms != "", .N],
    type = "percentage") %>%
    .[, sample_mean := round(n/N * 100, 1)] %>%
    .[, c("pop_mean", "pop_confint_low", "pop_confint_high") :=
        svyMean2(~symptom_respiratory=="yes", participant_data_design_ps %>% subset(symptoms != ""),
                 multiply = 100, digits=1, na.rm=T) %>%
        .[option == TRUE, c("mean", "ci95_low", "ci95_high")]] %>%
  .[, c("popu5_mean", "popu5_confint_low", "popu5_confint_high") :=
      svyMean2(~symptom_respiratory=="yes", participant_data_u5_design_ps %>% subset(symptoms != ""),
               multiply = 100, digits=1, na.rm=T) %>%
      .[option == TRUE, c("mean", "ci95_low", "ci95_high")]],
  data.table(variable_name = "Antibiotic use",
    n = participant_data_design_ps$variables[!is.na(antibiotics) & antibiotics == "yes", .N],
    N = participant_data_design_ps$variables[!is.na(antibiotics) & antibiotics != "", .N],
    type = "percentage") %>%
    .[, sample_mean := round(n/N * 100, 1)] %>%
    .[, c("pop_mean", "pop_confint_low", "pop_confint_high") :=
        svyMean2(~antibiotics=="yes", participant_data_design_ps %>% subset(antibiotics != ""),
                 multiply = 100, digits=1, na.rm=T) %>%
        .[option == TRUE, c("mean", "ci95_low", "ci95_high")]] %>%
    .[, c("popu5_mean", "popu5_confint_low", "popu5_confint_high") :=
        svyMean2(~antibiotics=="yes", participant_data_u5_design_ps%>% subset(antibiotics != ""),
                 multiply = 100, digits=1, na.rm=T) %>%
        .[option == TRUE, c("mean", "ci95_low", "ci95_high")]],
  data.table(variable_name = "Direct contacts",
    n = participant_data_design_ps$variables[!is.na(contacts_total_reported), .N],
    N = participant_data_design_ps$variables[!is.na(contacts_total_reported), .N],
    type = "percentage") %>%
    .[, sample_mean := round(participant_data_design_ps$variables[!is.na(contacts_total_reported), mean(contacts_total_reported)], 1)] %>%
    .[, c("pop_mean", "pop_confint_low", "pop_confint_high") :=
        svyMean2(~contacts_total_reported, participant_data_design_ps, digits=1, na.rm=T) %>%
        .[, c("mean", "ci95_low", "ci95_high")]] %>%
    .[, c("popu5_mean", "popu5_confint_low", "popu5_confint_high") :=
        svyMean2(~contacts_total_reported, participant_data_u5_design_ps, digits=1, na.rm=T) %>%
        .[, c("mean", "ci95_low", "ci95_high")]])
table1.2[, heading := "Potential risk factors"]

#############################################################################################################################
#' Section 3. Pneumococcal carriage prevalence
#############################################################################################################################
table1.3 = rbind(
  data.table(variable_name = "Pneumococcal carriage",
    n = participant_data_design_ps$variables[!is.na(pneu_carr_final) & pneu_carr_final == 1, .N],
    N = participant_data_design_ps$variables[!is.na(pneu_carr_final), .N],
    type = "percentage") %>%
    .[, sample_mean := round(n/N * 100, 1)] %>%
    .[, c("pop_mean", "pop_confint_low", "pop_confint_high") :=
        svyMean2(~pneu_carr_final, participant_data_design_ps, multiply = 100, digits=1, na.rm=T) %>%
        .[, c("mean", "ci95_low", "ci95_high")]] %>%
    .[, c("popu5_mean", "popu5_confint_low", "popu5_confint_high") :=
        svyMean2(~pneu_carr_final, participant_data_u5_design_ps, multiply = 100, digits=1, na.rm=T) %>%
        .[, c("mean", "ci95_low", "ci95_high")]],
  data.table(variable_name = "Non-encapsulated serotype carriers",
    n = participant_data_design_ps$variables[NT == 1, .N],
    N = participant_data_design_ps$variables[!is.na(NT) & pneu_carr_final == 1, .N],
    type = "percentage") %>%
    .[, sample_mean := round(n/N * 100, 1)] %>%
    .[, c("pop_mean", "pop_confint_low", "pop_confint_high") :=
        svyMean2(~NT, participant_data_design_ps %>% subset(pneu_carr_final == 1), multiply = 100, digits=1, na.rm=T) %>%
        .[, c("mean", "ci95_low", "ci95_high")]] %>%
    .[, c("popu5_mean", "popu5_confint_low", "popu5_confint_high") :=
        svyMean2(~NT, participant_data_u5_design_ps %>% subset(pneu_carr_final == 1), multiply = 100, digits=1, na.rm=T) %>%
        .[, c("mean", "ci95_low", "ci95_high")]])
table1.3[, heading := "Carriage prevalence"]  

#############################################################################################################################
#' Section 4. Proportion carriers with VT, using different definitions
#############################################################################################################################
table1.4 = rbind(
  data.table(variable_name = "PNEUMOSIL",
             n = participant_data_design_ps$variables[!is.na(pneumosil_carrier) & pneumosil_carrier == 1, .N],
             N = participant_data_design_ps$variables[!is.na(pneu_carr_final) & pneu_carr_final == 1, .N],
             type = "percentage") %>%
    .[, sample_mean := round(n/N * 100, 1)] %>%
    .[, c("pop_mean", "pop_confint_low", "pop_confint_high") :=
        svyMean2(~pneumosil_carrier, participant_data_design_ps %>% subset(!is.na(pneu_carr_final) & pneu_carr_final == 1),
                 multiply = 100, digits=1, na.rm=T) %>% .[, c("mean", "ci95_low", "ci95_high")]] %>%
    .[, c("popu5_mean", "popu5_confint_low", "popu5_confint_high") :=
        svyMean2(~pneumosil_carrier, participant_data_u5_design_ps %>% subset(!is.na(pneu_carr_final) & pneu_carr_final == 1),
                 multiply = 100, digits=1, na.rm=T) %>% .[, c("mean", "ci95_low", "ci95_high")]],
  data.table(variable_name = "Synflorix",
             n = participant_data_design_ps$variables[!is.na(pneumosil_carrier) & pcv10_carrier == 1, .N],
             N = participant_data_design_ps$variables[!is.na(pneu_carr_final) & pneu_carr_final == 1, .N],
             type = "percentage") %>%
    .[, sample_mean := round(n/N * 100, 1)] %>%
    .[, c("pop_mean", "pop_confint_low", "pop_confint_high") :=
        svyMean2(~pcv10_carrier, participant_data_design_ps %>% subset(!is.na(pneu_carr_final) & pneu_carr_final == 1),
                 multiply = 100, digits=1, na.rm=T) %>% .[, c("mean", "ci95_low", "ci95_high")]] %>%
    .[, c("popu5_mean", "popu5_confint_low", "popu5_confint_high") :=
        svyMean2(~pcv10_carrier, participant_data_u5_design_ps %>% subset(!is.na(pneu_carr_final) & pneu_carr_final == 1),
                 multiply = 100, digits=1, na.rm=T) %>% .[, c("mean", "ci95_low", "ci95_high")]],
  data.table(variable_name = "Prevenar 13",
             n = participant_data_design_ps$variables[!is.na(pneumosil_carrier) & pcv13_carrier == 1, .N],
             N = participant_data_design_ps$variables[!is.na(pneu_carr_final) & pneu_carr_final == 1, .N],
             type = "percentage") %>%
    .[, sample_mean := round(n/N * 100, 1)] %>%
    .[, c("pop_mean", "pop_confint_low", "pop_confint_high") :=
        svyMean2(~pcv13_carrier, participant_data_design_ps %>% subset(!is.na(pneu_carr_final) & pneu_carr_final == 1),
                 multiply = 100, digits=1, na.rm=T) %>% .[, c("mean", "ci95_low", "ci95_high")]] %>%
    .[, c("popu5_mean", "popu5_confint_low", "popu5_confint_high") :=
        svyMean2(~pcv13_carrier, participant_data_u5_design_ps %>% subset(!is.na(pneu_carr_final) & pneu_carr_final == 1),
                 multiply = 100, digits=1, na.rm=T) %>% .[, c("mean", "ci95_low", "ci95_high")]],
  data.table(variable_name = "Vaxneuvance",
             n = participant_data_design_ps$variables[!is.na(pneumosil_carrier) & pcv15_carrier == 1, .N],
             N = participant_data_design_ps$variables[!is.na(pneu_carr_final) & pneu_carr_final == 1, .N],
             type = "percentage") %>%
    .[, sample_mean := round(n/N * 100, 1)] %>%
    .[, c("pop_mean", "pop_confint_low", "pop_confint_high") :=
        svyMean2(~pcv15_carrier, participant_data_design_ps %>% subset(!is.na(pneu_carr_final) & pneu_carr_final == 1),
                 multiply = 100, digits=1, na.rm=T) %>% .[, c("mean", "ci95_low", "ci95_high")]] %>%
    .[, c("popu5_mean", "popu5_confint_low", "popu5_confint_high") :=
        svyMean2(~pcv15_carrier, participant_data_u5_design_ps %>% subset(!is.na(pneu_carr_final) & pneu_carr_final == 1),
                 multiply = 100, digits=1, na.rm=T) %>% .[, c("mean", "ci95_low", "ci95_high")]],
  data.table(variable_name = "Prevenar 20",
             n = participant_data_design_ps$variables[!is.na(pneumosil_carrier) & pcv20_carrier == 1, .N],
             N = participant_data_design_ps$variables[!is.na(pneu_carr_final) & pneu_carr_final == 1, .N],
             type = "percentage") %>%
    .[, sample_mean := round(n/N * 100, 1)] %>%
    .[, c("pop_mean", "pop_confint_low", "pop_confint_high") :=
        svyMean2(~pcv20_carrier, participant_data_design_ps %>% subset(!is.na(pneu_carr_final) & pneu_carr_final == 1),
                 multiply = 100, digits=1, na.rm=T) %>% .[, c("mean", "ci95_low", "ci95_high")]] %>%
    .[, c("popu5_mean", "popu5_confint_low", "popu5_confint_high") :=
        svyMean2(~pcv20_carrier, participant_data_u5_design_ps %>% subset(!is.na(pneu_carr_final) & pneu_carr_final == 1),
                 multiply = 100, digits=1, na.rm=T) %>% .[, c("mean", "ci95_low", "ci95_high")]])
table1.4[, heading := "Proportion of carriers with VT"] 

#############################################################################################################################
#' Section 5. Estimates for the final section of this table require a more detailed analyses, and are
#'            generated in 3_invasiveness_estimates.R
#############################################################################################################################

#############################################################################################################################
#' Section 6. Combine all sections in the table
#############################################################################################################################
table1_characteristics_and_prevalence =
  rbind(table1.1, table1.2, table1.3, table1.4, fill=T, use.names=T) %>%
  .[, c("heading", "variable_name", "n", "N", "sample_mean", "sample_confint_low", "sample_confint_high", "pop_mean",
        "pop_confint_low", "pop_confint_high", "popu5_mean", "popu5_confint_low", "popu5_confint_high", "type")]
table1_characteristics_and_prevalence[, n_format := ifelse(is.na(n), N, sprintf("%s/%s", n, N))]
table1_characteristics_and_prevalence[!is.na(sample_confint_low), sample_confint :=
                                          sprintf("%s - %s (IQR)", pmax(0, sample_confint_low), sample_confint_high)]
table1_characteristics_and_prevalence[, pop_confint :=
                                          sprintf("%s - %s", pmax(0, pop_confint_low), pop_confint_high)]
table1_characteristics_and_prevalence[, popu5_confint :=
                                          sprintf("%s - %s", pmax(0, popu5_confint_low), popu5_confint_high)]
table1_characteristics_and_prevalence[, sample_mean_format := ifelse(type == "percentage", sprintf("%s%%", sample_mean), sample_mean)]
table1_characteristics_and_prevalence[, pop_mean_format := ifelse(type == "percentage", sprintf("%s%%", pop_mean), pop_mean)]
table1_characteristics_and_prevalence[, popu5_mean_format := ifelse(type == "percentage", sprintf("%s%%", popu5_mean), popu5_mean)]
rm("table1.1", "table1.2", "table1.3", "table1.4")

#############################################################################################################################
#' Generate Table 1: "Sample characteristics and carriage prevalence."
#############################################################################################################################
if(make_table) table1_characteristics_and_prevalence[, c("heading", "variable_name", "n_format",
                                                           "sample_mean_format", "sample_confint",
                                                           "pop_mean_format", "pop_confint",
                                                           "popu5_mean_format", "popu5_confint")] %>%
  kblOut(booktabs=T, align=c("l","l","l","r","l","r","l","r","l"), linesep = "",
         col.names=c("", "", "total", "sample value", "", "population estimate", "", "population estimate (<5y)", ""),
         other_functions = list(function(x) collapse_rows(x, 1:2, row_group_label_position = "stack",
                                                          row_group_label_fonts = list(list(bold=T, italic=T),
                                                                                       list(bold=T, italic=T))),
                                function(x) kable_styling(x, latex_options = "scale_down")),
         out_name = "table1_characteristics_and_prevalence")
