#' some strata will be small - use these options for conservative estimates of the variance estimator
options(survey.lonely.psu="adjust")
options(survey.adjust.domain.lonely=TRUE)

#############################################################################################################################
#' Section 1. Assess differences
#' - No evidence of a difference, but sample sizes are relatively small
#############################################################################################################################  
#' Assess whether there is any difference in the odds of carriage between the two shipments
glm(pneu_carr_final~pilot, data=participant_data, family=stats::binomial()) %>%
  parseGLM("pilot") %>% .[]

#' Assess whether there is any difference in the number of carried serotypes, between two shipments
glm(num_sero_final~pilot, data=participant_data[pneu_carr_final == 1], family=stats::poisson()) %>%
  parseGLM("pilot") %>% .[]

#' Assess whether there is any difference in the log10-density of carried serotypes between the two shipments
glm(log10(pneu_dens_final)~pilot, data=subset(participant_data, pneu_dens_final>1), family=stats::gaussian) %>%
  parseGLM("pilot", exp = FALSE) %>% .[]

#############################################################################################################################
#' Section 2. Prevalence in pilot shipment
#############################################################################################################################
tableSB2.1 = rbind(
  data.table(variable_name = "Pneumococcal carriage",
             n = participant_data[pilot == 1 & !is.na(pneu_carr_final) & pneu_carr_final == 1, .N],
             #' Nb, we were not able to match the data from 2 of the swabs in the pilot shipment
             #' - this affects the population level estimates
             N = participant_data[pilot == 1 & !is.na(pneu_carr_final), .N], 
             type = "percentage") %>%
    .[, sample_mean := round(n/N * 100, 1)] %>%
    .[, c("pop_mean", "pop_confint_low", "pop_confint_high") :=
        svyMean2(~pneu_carr_final, participant_data_pilot_design_ps, multiply = 100, digits=1, na.rm=T) %>%
        .[, c("mean", "ci95_low", "ci95_high")]] %>%
    .[, c("popu5_mean", "popu5_confint_low", "popu5_confint_high") :=
        svyMean2(~pneu_carr_final, participant_data_u5_pilot_design_ps, multiply = 100, digits=1, na.rm=T) %>%
        .[, c("mean", "ci95_low", "ci95_high")]],
  {
    participant_data[, pneumosil_carriage := ifelse(is.na(pneu_carr_final), NA_integer_,
                              pneumosil_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))]
    participant_data_pilot_design_ps$variables[, pneumosil_carriage :=
                                           ifelse(is.na(pneu_carr_final), NA_integer_,
                                                  pneumosil_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))]
    participant_data_u5_pilot_design_ps$variables[, pneumosil_carriage :=
                                              ifelse(is.na(pneu_carr_final), NA_integer_,
                                                     pneumosil_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))]
    data.table(variable_name = "PNEUMOSIL-covered serotype carriers",
               n = participant_data[pilot == 1 & !is.na(pneumosil_type) & pneumosil_carriage == 1, .N],
               N = participant_data[pilot == 1 & !is.na(pneu_carr_final) & pneu_carr_final == 1, .N],
               type = "percentage") %>%
      .[, sample_mean := round(n/N * 100, 1)] %>%
      .[, c("pop_mean", "pop_confint_low", "pop_confint_high") :=
          svyMean2(~pneumosil_carriage, participant_data_pilot_design_ps %>% subset(pneu_carr_final == 1), multiply = 100, digits=1, na.rm=T) %>%
          .[1, c("mean", "ci95_low", "ci95_high")]] %>%
      .[, c("popu5_mean", "popu5_confint_low", "popu5_confint_high") :=
          svyMean2(~pneumosil_carriage, participant_data_u5_pilot_design_ps %>% subset(pneu_carr_final == 1), multiply = 100, digits=1, na.rm=T) %>%
          .[1, c("mean", "ci95_low", "ci95_high")]] %>% .[]},
  {
    participant_data[, pcv10_carriage :=
                       ifelse(is.na(pneu_carr_final), NA_integer_,
                              pcv10_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))]
    participant_data_pilot_design_ps$variables[, pcv10_carriage :=
                                           ifelse(is.na(pneu_carr_final), NA_integer_,
                                                  pcv10_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))];
    participant_data_u5_pilot_design_ps$variables[, pcv10_carriage :=
                                              ifelse(is.na(pneu_carr_final), NA_integer_,
                                                     pcv10_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))];
    
    data.table(variable_name = "Synflorix-covered serotype carriers",
               n = participant_data[pilot == 1 & !is.na(pcv10_carriage) & pcv10_carriage == 1, .N],
               N = participant_data[pilot == 1 & !is.na(pcv10_carriage) & pneu_carr_final == 1, .N],
               type = "percentage") %>%
      .[, sample_mean := round(n/N * 100, 1)] %>%
      .[, c("pop_mean", "pop_confint_low", "pop_confint_high") :=
          svyMean2(~pcv10_carriage, participant_data_pilot_design_ps %>% subset(pneu_carr_final == 1), multiply = 100, digits=1, na.rm=T) %>%
          .[1, c("mean", "ci95_low", "ci95_high")]] %>%
      .[, c("popu5_mean", "popu5_confint_low", "popu5_confint_high") :=
          svyMean2(~pcv10_carriage, participant_data_u5_pilot_design_ps %>% subset(pneu_carr_final == 1), multiply = 100, digits=1, na.rm=T) %>%
          .[1, c("mean", "ci95_low", "ci95_high")]] %>% .[]},
  {
    participant_data[, pcv13_carriage :=
                       ifelse(is.na(pneu_carr_final), NA_integer_,
                              pcv13_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))]
    participant_data_pilot_design_ps$variables[, pcv13_carriage :=
                                           ifelse(is.na(pneu_carr_final), NA_integer_,
                                                  pcv13_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))];
    participant_data_u5_pilot_design_ps$variables[, pcv13_carriage :=
                                              ifelse(is.na(pneu_carr_final), NA_integer_,
                                                     pcv13_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))];
    
    data.table(variable_name = "Prevenar 13-covered serotype carriers",
               n = participant_data[pilot == 1 & !is.na(pcv13_carriage) & pcv13_carriage == 1, .N],
               N = participant_data[pilot == 1 & !is.na(pcv13_carriage) & pneu_carr_final == 1, .N],
               type = "percentage") %>%
      .[, sample_mean := round(n/N * 100, 1)] %>%
      .[, c("pop_mean", "pop_confint_low", "pop_confint_high") :=
          svyMean2(~pcv13_carriage, participant_data_pilot_design_ps %>% subset(pneu_carr_final == 1), multiply = 100, digits=1, na.rm=T) %>%
          .[1, c("mean", "ci95_low", "ci95_high")]] %>%
      .[, c("popu5_mean", "popu5_confint_low", "popu5_confint_high") :=
          svyMean2(~pcv13_carriage, participant_data_u5_pilot_design_ps %>% subset(pneu_carr_final == 1), multiply = 100, digits=1, na.rm=T) %>%
          .[1, c("mean", "ci95_low", "ci95_high")]] %>% .[]},
  {
    participant_data[, pcv15_carriage :=
                       ifelse(is.na(pneu_carr_final), NA_integer_,
                              pcv15_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))]
    participant_data_pilot_design_ps$variables[, pcv15_carriage :=
                                           ifelse(is.na(pneu_carr_final), NA_integer_,
                                                  pcv15_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))];
    participant_data_u5_pilot_design_ps$variables[, pcv15_carriage :=
                                              ifelse(is.na(pneu_carr_final), NA_integer_,
                                                     pcv15_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))];
    
    data.table(variable_name = "Vaxneuvance-covered serotype carriers",
               n = participant_data[pilot == 1 & !is.na(pcv15_carriage) & pcv15_carriage == 1, .N],
               N = participant_data[pilot == 1 & !is.na(pcv15_carriage) & pneu_carr_final == 1, .N],
               type = "percentage") %>%
      .[, sample_mean := round(n/N * 100, 1)] %>%
      .[, c("pop_mean", "pop_confint_low", "pop_confint_high") :=
          svyMean2(~pcv15_carriage, participant_data_pilot_design_ps %>% subset(pneu_carr_final == 1), multiply = 100, digits=1, na.rm=T) %>%
          .[1, c("mean", "ci95_low", "ci95_high")]] %>%
      .[, c("popu5_mean", "popu5_confint_low", "popu5_confint_high") :=
          svyMean2(~pcv15_carriage, participant_data_u5_pilot_design_ps %>% subset(pneu_carr_final == 1), multiply = 100, digits=1, na.rm=T) %>%
          .[1, c("mean", "ci95_low", "ci95_high")]] %>% .[]},
  {
    participant_data[, pcv20_carriage :=
                       ifelse(is.na(pneu_carr_final), NA_integer_,
                              pcv20_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))]
    participant_data_pilot_design_ps$variables[, pcv20_carriage :=
                                           ifelse(is.na(pneu_carr_final), NA_integer_,
                                                  pcv20_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))];
    participant_data_u5_pilot_design_ps$variables[, pcv20_carriage :=
                                              ifelse(is.na(pneu_carr_final), NA_integer_,
                                                     pcv20_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))];
    
    data.table(variable_name = "Prevenar 20-covered serotype carriers",
               n = participant_data[pilot == 1 & !is.na(pcv20_carriage) & pcv20_carriage == 1, .N],
               N = participant_data[pilot == 1 & !is.na(pcv20_carriage) & pneu_carr_final == 1, .N],
               type = "percentage") %>%
      .[, sample_mean := round(n/N * 100, 1)] %>%
      .[, c("pop_mean", "pop_confint_low", "pop_confint_high") :=
          svyMean2(~pcv20_carriage, participant_data_pilot_design_ps %>% subset(pneu_carr_final == 1), multiply = 100, digits=1, na.rm=T) %>%
          .[1, c("mean", "ci95_low", "ci95_high")]] %>%
      .[, c("popu5_mean", "popu5_confint_low", "popu5_confint_high") :=
          svyMean2(~pcv20_carriage, participant_data_u5_pilot_design_ps %>% subset(pneu_carr_final == 1), multiply = 100, digits=1, na.rm=T) %>%
          .[1, c("mean", "ci95_low", "ci95_high")]] %>% .[]})
tableSB2.1[, heading := "Pilot shipment"]

#############################################################################################################################
#' Section 3. Prevalence in second shipment
#############################################################################################################################
tableSB2.2 = rbind(
  data.table(variable_name = "Pneumococcal carriage",
             n = participant_data[pilot == 0 & !is.na(pneu_carr_final) & pneu_carr_final == 1, .N],
             N = participant_data[pilot == 0 & !is.na(pneu_carr_final), .N], 
             type = "percentage") %>%
    .[, sample_mean := round(n/N * 100, 1)] %>%
    .[, c("pop_mean", "pop_confint_low", "pop_confint_high") :=
        svyMean2(~pneu_carr_final, participant_data_nopilot_design_ps, multiply = 100, digits=1, na.rm=T) %>%
        .[, c("mean", "ci95_low", "ci95_high")]] %>%
    .[, c("popu5_mean", "popu5_confint_low", "popu5_confint_high") :=
        svyMean2(~pneu_carr_final, participant_data_u5_nopilot_design_ps, multiply = 100, digits=1, na.rm=T) %>%
        .[, c("mean", "ci95_low", "ci95_high")]],
  {
    participant_data[, pneumosil_carriage :=
                       ifelse(is.na(pneu_carr_final), NA_integer_,
                              pneumosil_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))]
    participant_data_nopilot_design_ps$variables[, pneumosil_carriage :=
                                                 ifelse(is.na(pneu_carr_final), NA_integer_,
                                                        pneumosil_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))];
    participant_data_u5_nopilot_design_ps$variables[, pneumosil_carriage :=
                                                    ifelse(is.na(pneu_carr_final), NA_integer_,
                                                           pneumosil_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))];
    
    data.table(variable_name = "PNEUMOSIL-covered serotype carriers",
               n = participant_data[pilot == 0 & !is.na(pneumosil_type) & pneumosil_carriage == 1, .N],
               N = participant_data[pilot == 0 & !is.na(pneu_carr_final) & pneu_carr_final == 1, .N],
               type = "percentage") %>%
      .[, sample_mean := round(n/N * 100, 1)] %>%
      .[, c("pop_mean", "pop_confint_low", "pop_confint_high") :=
          svyMean2(~pneumosil_carriage, participant_data_nopilot_design_ps %>% subset(pneu_carr_final == 1), multiply = 100, digits=1, na.rm=T) %>%
          .[1, c("mean", "ci95_low", "ci95_high")]] %>%
      .[, c("popu5_mean", "popu5_confint_low", "popu5_confint_high") :=
          svyMean2(~pneumosil_carriage, participant_data_u5_nopilot_design_ps %>% subset(pneu_carr_final == 1), multiply = 100, digits=1, na.rm=T) %>%
          .[1, c("mean", "ci95_low", "ci95_high")]] %>% .[]},
  {
    participant_data[, pcv10_carriage :=
                       ifelse(is.na(pneu_carr_final), NA_integer_,
                              pcv10_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))]
    participant_data_nopilot_design_ps$variables[, pcv10_carriage :=
                                                 ifelse(is.na(pneu_carr_final), NA_integer_,
                                                        pcv10_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))];
    participant_data_u5_nopilot_design_ps$variables[, pcv10_carriage :=
                                                    ifelse(is.na(pneu_carr_final), NA_integer_,
                                                           pcv10_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))];
    
    data.table(variable_name = "Synflorix-covered serotype carriers",
               n = participant_data[pilot == 0 & !is.na(pcv10_carriage) & pcv10_carriage == 1, .N],
               N = participant_data[pilot == 0 & !is.na(pcv10_carriage) & pneu_carr_final == 1, .N],
               type = "percentage") %>%
      .[, sample_mean := round(n/N * 100, 1)] %>%
      .[, c("pop_mean", "pop_confint_low", "pop_confint_high") :=
          svyMean2(~pcv10_carriage, participant_data_nopilot_design_ps %>% subset(pneu_carr_final == 1), multiply = 100, digits=1, na.rm=T) %>%
          .[1, c("mean", "ci95_low", "ci95_high")]] %>%
      .[, c("popu5_mean", "popu5_confint_low", "popu5_confint_high") :=
          svyMean2(~pcv10_carriage, participant_data_u5_nopilot_design_ps %>% subset(pneu_carr_final == 1), multiply = 100, digits=1, na.rm=T) %>%
          .[1, c("mean", "ci95_low", "ci95_high")]] %>% .[]},
  {
    participant_data[, pcv13_carriage :=
                       ifelse(is.na(pneu_carr_final), NA_integer_,
                              pcv13_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))]
    participant_data_nopilot_design_ps$variables[, pcv13_carriage :=
                                                 ifelse(is.na(pneu_carr_final), NA_integer_,
                                                        pcv13_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))];
    participant_data_u5_nopilot_design_ps$variables[, pcv13_carriage :=
                                                    ifelse(is.na(pneu_carr_final), NA_integer_,
                                                           pcv13_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))];
    
    data.table(variable_name = "Prevenar 13-covered serotype carriers",
               n = participant_data[pilot == 0 & !is.na(pcv13_carriage) & pcv13_carriage == 1, .N],
               N = participant_data[pilot == 0 & !is.na(pcv13_carriage) & pneu_carr_final == 1, .N],
               type = "percentage") %>%
      .[, sample_mean := round(n/N * 100, 1)] %>%
      .[, c("pop_mean", "pop_confint_low", "pop_confint_high") :=
          svyMean2(~pcv13_carriage, participant_data_nopilot_design_ps %>% subset(pneu_carr_final == 1), multiply = 100, digits=1, na.rm=T) %>%
          .[1, c("mean", "ci95_low", "ci95_high")]] %>%
      .[, c("popu5_mean", "popu5_confint_low", "popu5_confint_high") :=
          svyMean2(~pcv13_carriage, participant_data_u5_nopilot_design_ps %>% subset(pneu_carr_final == 1), multiply = 100, digits=1, na.rm=T) %>%
          .[1, c("mean", "ci95_low", "ci95_high")]] %>% .[]},
  {
    participant_data[, pcv15_carriage :=
                       ifelse(is.na(pneu_carr_final), NA_integer_,
                              pcv15_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))]
    participant_data_nopilot_design_ps$variables[, pcv15_carriage :=
                                                 ifelse(is.na(pneu_carr_final), NA_integer_,
                                                        pcv15_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))];
    participant_data_u5_nopilot_design_ps$variables[, pcv15_carriage :=
                                                    ifelse(is.na(pneu_carr_final), NA_integer_,
                                                           pcv15_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))];
    
    data.table(variable_name = "Vaxneuvance-covered serotype carriers",
               n = participant_data[pilot == 0 & !is.na(pcv15_carriage) & pcv15_carriage == 1, .N],
               N = participant_data[pilot == 0 & !is.na(pcv15_carriage) & pneu_carr_final == 1, .N],
               type = "percentage") %>%
      .[, sample_mean := round(n/N * 100, 1)] %>%
      .[, c("pop_mean", "pop_confint_low", "pop_confint_high") :=
          svyMean2(~pcv15_carriage, participant_data_nopilot_design_ps %>% subset(pneu_carr_final == 1), multiply = 100, digits=1, na.rm=T) %>%
          .[1, c("mean", "ci95_low", "ci95_high")]] %>%
      .[, c("popu5_mean", "popu5_confint_low", "popu5_confint_high") :=
          svyMean2(~pcv15_carriage, participant_data_u5_nopilot_design_ps %>% subset(pneu_carr_final == 1), multiply = 100, digits=1, na.rm=T) %>%
          .[1, c("mean", "ci95_low", "ci95_high")]] %>% .[]},
  {
    participant_data[, pcv20_carriage :=
                       ifelse(is.na(pneu_carr_final), NA_integer_,
                              pcv20_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))]
    participant_data_nopilot_design_ps$variables[, pcv20_carriage :=
                                                 ifelse(is.na(pneu_carr_final), NA_integer_,
                                                        pcv20_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))];
    participant_data_u5_nopilot_design_ps$variables[, pcv20_carriage :=
                                                    ifelse(is.na(pneu_carr_final), NA_integer_,
                                                           pcv20_type %in% c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT"))];
    
    data.table(variable_name = "Prevenar 20-covered serotype carriers",
               n = participant_data[pilot == 0 & !is.na(pcv20_carriage) & pcv20_carriage == 1, .N],
               N = participant_data[pilot == 0 & !is.na(pcv20_carriage) & pneu_carr_final == 1, .N],
               type = "percentage") %>%
      .[, sample_mean := round(n/N * 100, 1)] %>%
      .[, c("pop_mean", "pop_confint_low", "pop_confint_high") :=
          svyMean2(~pcv20_carriage, participant_data_nopilot_design_ps %>% subset(pneu_carr_final == 1), multiply = 100, digits=1, na.rm=T) %>%
          .[1, c("mean", "ci95_low", "ci95_high")]] %>%
      .[, c("popu5_mean", "popu5_confint_low", "popu5_confint_high") :=
          svyMean2(~pcv20_carriage, participant_data_u5_nopilot_design_ps %>% subset(pneu_carr_final == 1), multiply = 100, digits=1, na.rm=T) %>%
          .[1, c("mean", "ci95_low", "ci95_high")]] %>% .[]})
tableSB2.2[, heading := "Second shipment"]

#############################################################################################################################
#' Section 4. Generate Supplemental Table B2. Carriage prevalence in pilot and second shipment.
#############################################################################################################################
tableSB2_shipments_prevalence =
  rbind(tableSB2.1, tableSB2.2, fill=T, use.names=T) %>%
  .[, c("heading", "variable_name", "n", "N", "sample_mean", "pop_mean", "pop_confint_low", "pop_confint_high",
        "popu5_mean", "popu5_confint_low", "popu5_confint_high", "type")]
tableSB2_shipments_prevalence[, n_format := ifelse(is.na(n), N, sprintf("%s/%s", n, N))]
tableSB2_shipments_prevalence[, pop_confint :=
                                          sprintf("%s - %s", pmax(0, pop_confint_low), pop_confint_high)]
tableSB2_shipments_prevalence[, popu5_confint :=
                                          sprintf("%s - %s", pmax(0, popu5_confint_low), popu5_confint_high)]
tableSB2_shipments_prevalence[, sample_mean_format := ifelse(type == "percentage", sprintf("%s%%", sample_mean), sample_mean)]
tableSB2_shipments_prevalence[, pop_mean_format := ifelse(type == "percentage", sprintf("%s%%", pop_mean), pop_mean)]
tableSB2_shipments_prevalence[, popu5_mean_format := ifelse(type == "percentage", sprintf("%s%%", popu5_mean), popu5_mean)]
rm("tableSB2.1", "tableSB2.2")

#' Write the table to the output folder
if(make_table) tableSB2_shipments_prevalence[, c("heading", "variable_name", "n_format",
                                                 "sample_mean_format", "pop_mean_format", "pop_confint",
                                                 "popu5_mean_format", "popu5_confint")] %>%
  kblOut(booktabs=T, align=c("l","l","l","r","r","l","r","l"), linesep = "",
         col.names=c("", "", "total", "sample value", "population estimate", "", "population estimate (<5y)", ""),
         other_functions = list(function(x) collapse_rows(x, 1:2, row_group_label_position = "stack",
                                                          row_group_label_fonts = list(list(bold=T, italic=T),
                                                                                       list(bold=T, italic=T))),
                                function(x) kable_styling(x, latex_options = "scale_down")),
         out_name = "tableSB2_shipments_prevalence")
