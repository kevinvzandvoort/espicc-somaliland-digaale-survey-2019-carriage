#' Some values in this script are read from other files. These values are used to calculate the post-stratification weights.
#' They are equivalent to those in the https://github.com/kevinvzandvoort/espicc-somaliland-digaale-survey-2019 repository.
#' Download the kevinvzandvoort/espicc-somaliland-digaale-survey-2019 repository to see calculations for these values.

#######################################################
#' Section 1. Process swab data
#######################################################
#' 1. generate variables (pcv20, pcv20_15, pcv20_15_13, pcv20_15_13_10, pcv20_15_13_pneumosil, pcv20_15_13_10_pneumosil, vaccine)
#'    to assess whether each serotype was included PCV20 (Prevenar 20), PCV15 (Vaxneuvance), PCV13 (Prevenar 13),
#'    PCV10 (Synflorix), and/or Pneumosil
#' 2. use these variables to assess whether each individual/swab carried serotypes of each specific vaccine
participant_carriage_summary = lab_data %>%
  melt(measure.vars = list(serotype = c("sero1_final", "sero2_final", "sero3_final", "sero4_final", "sero5_final"),
                           relative_abundance = c("ra_sero1_final", "ra_sero2_final", "ra_sero3_final", "ra_sero4_final", "ra_sero5_final"))) %>%
  .[, c("pid", "serotype", "variable", "pneu_carr_final", "pneu_dens_final", "num_sero_final")] %>%
  .[, pcv20 := serotype %in% serotypes_by_vaccine[pcv20==1 & pcv15==0 & pcv13==0 & pcv10==0 & pneumosil==0, serotype]] %>%
  .[, pcv20_15 := serotype %in% serotypes_by_vaccine[pcv20==1 & pcv15==1 & pcv13==0 & pcv10==0 & pneumosil==0, serotype]] %>%
  .[, pcv20_15_13 := serotype %in% serotypes_by_vaccine[pcv20==1 & pcv15==1 & pcv13==1 & pcv10==0 & pneumosil==0, serotype]] %>%
  .[, pcv20_15_13_10 := serotype %in% serotypes_by_vaccine[pcv20==1 & pcv15==1 & pcv13==1 & pcv10==1 & pneumosil==0, serotype]] %>%
  .[, pcv20_15_13_pneumosil := serotype %in% serotypes_by_vaccine[pcv20==1 & pcv15==1 & pcv13==1 & pcv10==0 & pneumosil==1, serotype]] %>%
  .[, pcv20_15_13_10_pneumosil := serotype %in% serotypes_by_vaccine[pcv20==1 & pcv15==1 & pcv13==1 & pcv10==1 & pneumosil==1, serotype]] %>%
  .[, vaccine := ifelse(pcv20_15_13_10_pneumosil, "pcv20_15_13_10_pneumosil",
                        ifelse(pcv20_15_13_10, "pcv20_15_13_10",
                               ifelse(pcv20_15_13_pneumosil, "pcv20_15_13_pneumosil",
                                      ifelse(pcv20_15_13, "pcv20_15_13",
                                             ifelse(pcv20_15, "pcv20_15",
                                                    ifelse(pcv20, "pcv20",
                                                           ifelse(!is.na(serotype) & serotype != "", "NVT", NA_character_)))))))] %>%
  .[, vaccine := ifelse(grepl("NT", serotype, fixed=TRUE), "NT", vaccine)] %>%
  dcast(num_sero_final+pneu_carr_final+pneu_dens_final+pid~vaccine) %>%
  .[, pneumosil_type := ifelse(NT>0 & (NVT + pcv20_15_13_10 + pcv20_15_13 + pcv20_15 + pcv20)>0 & (pcv20_15_13_pneumosil + pcv20_15_13_10_pneumosil) > 0, "VT + NVT + NT",
                               ifelse(NT==0 & (NVT + pcv20_15_13_10 + pcv20_15_13 + pcv20_15 + pcv20)>0 & (pcv20_15_13_pneumosil + pcv20_15_13_10_pneumosil) > 0, "VT + NVT",
                                      ifelse(NT>0 & (NVT + pcv20_15_13_10 + pcv20_15_13 + pcv20_15 + pcv20)==0 & (pcv20_15_13_pneumosil + pcv20_15_13_10_pneumosil) > 0, "VT + NT",
                                             ifelse(NT==0 & (NVT + pcv20_15_13_10 + pcv20_15_13 + pcv20_15 + pcv20)==0 & (pcv20_15_13_pneumosil + pcv20_15_13_10_pneumosil) > 0, "VT",
                                                    ifelse(NT>0 & (NVT + pcv20_15_13_10 + pcv20_15_13 + pcv20_15 + pcv20)>0 & (pcv20_15_13_pneumosil + pcv20_15_13_10_pneumosil) == 0, "NVT + NT",
                                                           ifelse(NT==0 & (NVT + pcv20_15_13_10 + pcv20_15_13 + pcv20_15 + pcv20)>0 & (pcv20_15_13_pneumosil + pcv20_15_13_10_pneumosil) == 0, "NVT",
                                                                  ifelse(NT>0 & (NVT + pcv20_15_13_10 + pcv20_15_13 + pcv20_15 + pcv20)==0 & (pcv20_15_13_pneumosil + pcv20_15_13_10_pneumosil) == 0, "NT", "S")))))))] %>%
  .[, pcv10_type := ifelse(NT>0 & (NVT + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15 + pcv20)>0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil) > 0, "VT + NVT + NT",
                           ifelse(NT==0 & (NVT + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15 + pcv20)>0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil) > 0, "VT + NVT",
                                  ifelse(NT>0 & (NVT + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15 + pcv20)==0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil) > 0, "VT + NT",
                                         ifelse(NT==0 & (NVT + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15 + pcv20)==0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil) > 0, "VT",
                                                ifelse(NT>0 & (NVT + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15 + pcv20)>0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil) == 0, "NVT + NT",
                                                       ifelse(NT==0 & (NVT + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15 + pcv20)>0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil) == 0, "NVT",
                                                              ifelse(NT>0 & (NVT + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15 + pcv20)==0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil) == 0, "NT", "S")))))))] %>%
  .[, pcv13_type := ifelse(NT>0 & (NVT + pcv20_15 + pcv20)>0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13) > 0, "VT + NVT + NT",
                           ifelse(NT==0 & (NVT + pcv20_15 + pcv20)>0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13) > 0, "VT + NVT",
                                  ifelse(NT>0 & (NVT + pcv20_15 + pcv20)==0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13) > 0, "VT + NT",
                                         ifelse(NT==0 & (NVT + pcv20_15 + pcv20)==0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13) > 0, "VT",
                                                ifelse(NT>0 & (NVT + pcv20_15 + pcv20)>0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13) == 0, "NVT + NT",
                                                       ifelse(NT==0 & (NVT + pcv20_15 + pcv20)>0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13) == 0, "NVT",
                                                              ifelse(NT>0 & (NVT + pcv20_15 + pcv20)==0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13) == 0, "NT", "S")))))))] %>%
  .[, pcv15_type := ifelse(NT>0 & (NVT + pcv20)>0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15) > 0, "VT + NVT + NT",
                           ifelse(NT==0 & (NVT + pcv20)>0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15) > 0, "VT + NVT",
                                  ifelse(NT>0 & (NVT + pcv20)==0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15) > 0, "VT + NT",
                                         ifelse(NT==0 & (NVT + pcv20)==0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15) > 0, "VT",
                                                ifelse(NT>0 & (NVT + pcv20)>0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15) == 0, "NVT + NT",
                                                       ifelse(NT==0 & (NVT + pcv20)>0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15) == 0, "NVT",
                                                              ifelse(NT>0 & (NVT + pcv20)==0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15) == 0, "NT", "S")))))))] %>%
  .[, pcv20_type := ifelse(NT>0 & (NVT)>0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15 + pcv20) > 0, "VT + NVT + NT",
                           ifelse(NT==0 & (NVT)>0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15 + pcv20) > 0, "VT + NVT",
                                  ifelse(NT>0 & (NVT)==0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15 + pcv20) > 0, "VT + NT",
                                         ifelse(NT==0 & (NVT)==0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15 + pcv20) > 0, "VT",
                                                ifelse(NT>0 & (NVT)>0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15 + pcv20) == 0, "NVT + NT",
                                                       ifelse(NT==0 & (NVT)>0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15 + pcv20) == 0, "NVT",
                                                              ifelse(NT>0 & (NVT)==0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15 + pcv20) == 0, "NT", "S")))))))] %>%
  .[is.na(pneu_carr_final), c("NT") := .(NA_integer_)]

participant_carriage_summary[, pneumosil_carrier := ifelse(is.na(pneu_carr_final), NA_integer_, pneumosil_type %in% c("VT + NVT + NT", "VT + NVT", "VT + NT", "VT"))]
participant_carriage_summary[, pcv10_carrier := ifelse(is.na(pneu_carr_final), NA_integer_, pcv10_type %in% c("VT + NVT + NT", "VT + NVT", "VT + NT", "VT"))]
participant_carriage_summary[, pcv13_carrier := ifelse(is.na(pneu_carr_final), NA_integer_, pcv13_type %in% c("VT + NVT + NT", "VT + NVT", "VT + NT", "VT"))]
participant_carriage_summary[, pcv15_carrier := ifelse(is.na(pneu_carr_final), NA_integer_, pcv15_type %in% c("VT + NVT + NT", "VT + NVT", "VT + NT", "VT"))]
participant_carriage_summary[, pcv20_carrier := ifelse(is.na(pneu_carr_final), NA_integer_, pcv20_type %in% c("VT + NVT + NT", "VT + NVT", "VT + NT", "VT"))]

#' Same variable assignment as above, but only for the dominant serotype in a swab
participant_carriage_dominant_summary = lab_data %>%
  .[, c("pid", "sero1_final", "pneu_carr_final")] %>%
  .[, pcv20 := sero1_final %in% serotypes_by_vaccine[pcv20==1 & pcv15==0 & pcv13==0 & pcv10==0 & pneumosil==0, serotype]] %>%
  .[, pcv20_15 := sero1_final %in% serotypes_by_vaccine[pcv20==1 & pcv15==1 & pcv13==0 & pcv10==0 & pneumosil==0, serotype]] %>%
  .[, pcv20_15_13 := sero1_final %in% serotypes_by_vaccine[pcv20==1 & pcv15==1 & pcv13==1 & pcv10==0 & pneumosil==0, serotype]] %>%
  .[, pcv20_15_13_10 := sero1_final %in% serotypes_by_vaccine[pcv20==1 & pcv15==1 & pcv13==1 & pcv10==1 & pneumosil==0, serotype]] %>%
  .[, pcv20_15_13_pneumosil := sero1_final %in% serotypes_by_vaccine[pcv20==1 & pcv15==1 & pcv13==1 & pcv10==0 & pneumosil==1, serotype]] %>%
  .[, pcv20_15_13_10_pneumosil := sero1_final %in% serotypes_by_vaccine[pcv20==1 & pcv15==1 & pcv13==1 & pcv10==1 & pneumosil==1, serotype]] %>%
  .[, vaccine := ifelse(pcv20_15_13_10_pneumosil, "pcv20_15_13_10_pneumosil",
                        ifelse(pcv20_15_13_pneumosil, "pcv20_15_13_pneumosil",
                               ifelse(pcv20_15_13_pneumosil, "pcv20_15_13_pneumosil",
                                      ifelse(pcv20_15_13_10, "pcv20_15_13_10",
                                             ifelse(pcv20_15_13, "pcv20_15_13",
                                                    ifelse(pcv20_15, "pcv20_15",
                                                           ifelse(pcv20, "pcv20",
                                                                  ifelse(!is.na(pneu_carr_final) & pneu_carr_final == 1, "NVT",
                                                                         ifelse(!is.na(pneu_carr_final) & pneu_carr_final == 0, "S", NA_character_)))))))))] %>%
  .[, vaccine := ifelse(grepl("NT", sero1_final, fixed=TRUE), "NT", vaccine)] %>%
  dcast(pneu_carr_final+pid~vaccine, fun.aggregate=length) %>%
  #Add as no records are present
  .[, pcv20_15 := 0] %>%
  .[, dom_pneumosil_type := ifelse(NT==0 & (NVT + pcv20_15_13_10 + pcv20_15_13 + pcv20_15 + pcv20)==0 & (pcv20_15_13_pneumosil + pcv20_15_13_10_pneumosil) > 0, "VT",
                                   ifelse(NT==0 & (NVT + pcv20_15_13_10 + pcv20_15_13 + pcv20_15 + pcv20)>0 & (pcv20_15_13_pneumosil + pcv20_15_13_10_pneumosil) == 0, "NVT",
                                          ifelse(NT>0 & (NVT + pcv20_15_13_10 + pcv20_15_13 + pcv20_15 + pcv20)==0 & (pcv20_15_13_pneumosil + pcv20_15_13_10_pneumosil) == 0, "NT", "S")))] %>%
  .[, dom_pcv10_type := ifelse(NT==0 & (NVT + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15 + pcv20)==0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil) > 0, "VT",
                               ifelse(NT==0 & (NVT + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15 + pcv20)>0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil) == 0, "NVT",
                                      ifelse(NT>0 & (NVT + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15 + pcv20)==0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil) == 0, "NT", "S")))] %>%
  .[, dom_pcv13_type := ifelse(NT==0 & (NVT + pcv20_15 + pcv20)==0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13) > 0, "VT",
                               ifelse(NT==0 & (NVT + pcv20_15 + pcv20)>0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13) == 0, "NVT",
                                      ifelse(NT>0 & (NVT + pcv20_15 + pcv20)==0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13) == 0, "NT", "S")))] %>%
  .[, dom_pcv15_type := ifelse(NT==0 & (NVT + pcv20)==0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15) > 0, "VT",
                               ifelse(NT==0 & (NVT + pcv20)>0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15) == 0, "NVT",
                                      ifelse(NT>0 & (NVT + pcv20)==0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15) == 0, "NT", "S")))] %>%
  .[, dom_pcv20_type := ifelse(NT==0 & (NVT)==0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15 + pcv20) > 0, "VT",
                               ifelse(NT==0 & (NVT)>0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15 + pcv20) == 0, "NVT",
                                      ifelse(NT>0 & (NVT)==0 & (pcv20_15_13_10 + pcv20_15_13_10_pneumosil + pcv20_15_13_pneumosil + pcv20_15_13 + pcv20_15 + pcv20) == 0, "NT", "S")))] %>%
  .[, dom_NT := NT] %>% .[, -"NT"] %>%
  .[is.na(pneu_carr_final), c("dom_pneumosil_type", "dom_pcv10_type", "dom_pcv13_type", "dom_pcv15_type", "dom_pcv20_type", "dom_NT") := .(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_integer_)] %>%
  .[, c("pid", "dom_pneumosil_type", "dom_pcv10_type", "dom_pcv13_type", "dom_pcv15_type", "dom_pcv20_type", "dom_NT")]

#' add these new variables to the other participant data
participant_data = participant_data %>%
  merge(participant_carriage_summary, by="pid", all.x=TRUE) %>%
  merge(participant_carriage_dominant_summary, by="pid", all.x=TRUE)

#######################################################
#' Section 2. Process symptom data
#######################################################
participant_data[, symptom_headache := ifelse(is.na(symptoms) | symptoms == "", NA_integer_, grepl("Headache", symptoms))]
participant_data[, symptom_diarrhoea := ifelse(is.na(symptoms) | symptoms == "", NA_integer_, grepl("Diarrhoea", symptoms))]
participant_data[, symptom_fever_nochills := ifelse(is.na(symptoms) | symptoms == "", NA_integer_, grepl("Fever without chills", symptoms))]
participant_data[, symptom_fever_chills := ifelse(is.na(symptoms) | symptoms == "", NA_integer_, grepl("Fever with chills", symptoms))]
participant_data[, symptom_fever := ifelse(is.na(symptoms) | symptoms == "", NA_integer_, grepl("Fever", symptoms))]
participant_data[, symptom_cough := ifelse(is.na(symptoms) | symptoms == "", NA_integer_, grepl("Cough", symptoms))]
participant_data[, symptom_runnynose := ifelse(is.na(symptoms) | symptoms == "", NA_integer_, grepl("Runny nose", symptoms))]
participant_data[, symptom_sneezing := ifelse(is.na(symptoms) | symptoms == "", NA_integer_, grepl("Sneezing", symptoms))]
participant_data[, symptom_breathingdifficulty := ifelse(is.na(symptoms) | symptoms == "", NA_integer_, grepl("Difficulty breating", symptoms))]
participant_data[, symptom_sorethroat := ifelse(is.na(symptoms) | symptoms == "", NA_integer_, grepl("Sore throat", symptoms))]
participant_data[, symptom_wheezing := ifelse(is.na(symptoms) | symptoms == "", NA_integer_, grepl("Wheezing", symptoms))]
participant_data[, symptom_respiratory := (symptom_fever + symptom_cough + symptom_runnynose + symptom_sneezing + symptom_breathingdifficulty + symptom_sorethroat + symptom_wheezing) >= 1]
participant_data[, symptom_unknown := ifelse(is.na(symptoms) | symptoms == "", NA_integer_, grepl("Don't know", symptoms))]
participant_data[, symptom_none := ifelse(is.na(symptoms) | symptoms == "", NA_integer_, grepl("None", symptoms))]
participant_data[, symptom_other := ifelse(is.na(symptoms) | symptoms == "", NA_integer_, grepl("other", symptoms))]
participant_data[symptom_none == TRUE & symptom_cough == TRUE, symptom_none := FALSE]

#######################################################
#' Section 3. Process participant substance use data
#######################################################
participant_data[, substance_tobacco := ifelse(is.na(substance_use) | substance_use == "", NA_integer_, (grepl("Tobacco", substance_use) | grepl("Shisha", substance_use)))]
participant_data[, substance_khat := ifelse(is.na(substance_use) | substance_use == "", NA_integer_, grepl("Tobacco", substance_use))]
participant_data[, abs := ifelse(abs == "", NA_integer_, abs == "yes")]

#######################################################
#' Section 4. Process comorbidity data
#######################################################
participant_data[, morbidity_pneumonia := ifelse(is.na(morbidities) | morbidities == "", NA_integer_, grepl("Pneumonia", morbidities))]
participant_data[, morbidity_copd := ifelse(is.na(morbidities) | morbidities == "", NA_integer_, grepl("COPD", morbidities))]
participant_data[, morbidity_sicklecell := ifelse(is.na(morbidities) | morbidities == "", NA_integer_, grepl("Sickle", morbidities))]
participant_data[, morbidity_asthma := ifelse(is.na(morbidities) | morbidities == "", NA_integer_, grepl("Asthma", morbidities))]
participant_data[, morbidity_diabetes := ifelse(is.na(morbidities) | morbidities == "", NA_integer_, grepl("Diabetes", morbidities))]
participant_data[, pneumonia_6m := morbidity_pneumonia]
participant_data[pneumonia_6m == TRUE, pneumonia_6m := morbidities_pneumonia_time == "Less than six months ago"]

#######################################################
#' Section 5. Process shelter quality data
#######################################################
participant_data[, house_draft := factor(house_draft, c("no", "yes"))]
participant_data[, house_leakage := factor(house_leakage, c("no", "yes"))]
participant_data[, house_ventilation := factor(house_ventilation, c("no", "yes", "cook outside"))]

#######################################################
#' Section 6. Process age groups
#######################################################
participant_data[is.na(participant_age_month), participant_age_month := participant_age_y * 12]
participant_data[, participant_age_days := round(participant_age_month * (365/12))]

#' Age groups that were used in the sampling
age_groups_sample <- data.table(
  age_from = c(0, 2, 6, 15, 30, 50),
  age_to = c(1, 5, 14, 29, 49, 100)) %>%
  .[, name := paste0(age_from, "-", age_to)] %>%
  .[age_to == 1, name := "<2"] %>%
  .[age_to == 100, name := "50+"] %>% .[]

for(i in 1:nrow(age_groups_sample)){
  participant_data[participant_age_y >= age_groups_sample[i, age_from] & participant_age_y <= age_groups_sample[i, age_to],
                   agegrp := factor(
                     age_groups_sample[i, name], levels=age_groups_sample[, name])]
  participant_data[participant_age_y >= age_groups_sample[i, age_from] & participant_age_y <= age_groups_sample[i, age_to],
                   participant_age_group_sample := factor(
                     age_groups_sample[i, name], levels=age_groups_sample[, name])]
}

#######################################################
#' Section 7. Process contact data
#######################################################
#' Indirect contacts (contacts_other_est) are estimated as the median of the range in the selected category
participant_data = participant_data %>%
  .[contacts_other_est == "0", contacts_other_median := 0] %>%
  .[contacts_other_est == "1-2", contacts_other_median := median(1,2)] %>%
  .[contacts_other_est == "3-5", contacts_other_median := median(3,5)] %>%
  .[contacts_other_est == "6-10", contacts_other_median := median(6,10)] %>%
  .[contacts_other_est == "11-20", contacts_other_median := median(11,20)] %>%
  .[contacts_other_est == ">20", contacts_other_median := 20]
colnames(participant_data)[which(colnames(participant_data) == "contacts_total_recorded")] = "contacts_total_reported"
colnames(participant_data)[which(colnames(participant_data) == "contacts_physical_total_recorded")] = "contacts_physical_total_reported"

#######################################################
#' Section 8. Process nutrition data
#######################################################
#' zscorer package requires sex to be coded as 1 for males, and 2 for females
participant_data[, participant_sex_var := as.numeric(factor(participant_sex, levels=c("male", "female")))]
#' zscorer package requires standing to be coded as 1 for yes (measuring height), and 2 for no (measuring length)
participant_data[, standing := as.numeric(factor(height_type, levels=c("height", "length")))]
participant_data[is.na(height_type), standing := 3]

#' Use zscorer to estimate wfa, hfa, wfh, and mfa z-scores
wfa = zscorer::addWGSR(participant_data[, c("pid", "participant_sex_var", "weight", "participant_age_days", "standing")],
                       sex = "participant_sex_var", firstPart = "weight", secondPart = "participant_age_days",
                       index = "wfa", standing = "standing")
participant_data = merge(participant_data, wfa[,c("pid","wfaz")], by="pid")

hfa = zscorer::addWGSR(participant_data[, c("pid", "participant_sex_var", "height", "participant_age_days", "standing")],
                       sex = "participant_sex_var", firstPart = "height", secondPart = "participant_age_days",
                       index = "hfa", standing = "standing")
participant_data = merge(participant_data, hfa[, c("pid","hfaz")], by="pid")

wfh = zscorer::addWGSR(participant_data[, c("pid", "participant_sex_var", "height", "weight", "standing")],
                       sex = "participant_sex_var", firstPart = "weight", secondPart = "height", index = "wfh",
                       standing = "standing")
participant_data = merge(participant_data, wfh[,c("pid","wfhz")], by="pid")

mfa = zscorer::addWGSR(participant_data[, c("pid", "participant_sex_var", "muac", "participant_age_days", "standing")],
                       sex = "participant_sex_var", firstPart = "muac", secondPart = "participant_age_days",
                       index = "mfa", standing = "standing")
participant_data = merge(participant_data, mfa[,c("pid","mfaz")], by="pid")
rm("wfa", "hfa", "wfh", "mfa")

#' recode z-scores as categories
participant_data[wfaz > -2, weight_for_age := 1]
participant_data[wfaz <= -2, weight_for_age := 2]
participant_data[wfaz <= -3, weight_for_age := 3]
participant_data[, weight_for_age := factor(weight_for_age, 1:3, c("Not underweight (z > -2)", "Underweight (z <= -2)",
                                                                   "Underweight (z <= -3)"))]
participant_data[hfaz > -2, height_for_age := 1]
participant_data[hfaz <= -2, height_for_age := 2]
participant_data[hfaz <= -3, height_for_age := 3]
participant_data[, height_for_age := factor(height_for_age, 1:3, c("Not stunted (z > -2)", "Stunted (z <= -2)",
                                                                   "Severely stunted (z <= -3)"))]
participant_data[wfhz > -2, weight_for_height := 1]
participant_data[wfhz <= -2, weight_for_height := 2]
participant_data[wfhz <= -3, weight_for_height := 3]
participant_data[, weight_for_height := factor(weight_for_height, 1:3, c("Not wasted (z > -2)", "Wasted (z <= -2)",
                                                                         "Severely wasted (z <= -3)"))]
participant_data[(muac * 10) >= 125, muac_level := 1]
participant_data[(muac * 10) < 125, muac_level := 2]
participant_data[(muac * 10) < 115, muac_level := 3]
participant_data[, muac_level := factor(muac_level, 1:3, c("Not wasted (>= 125mm)", "Wasted (< -125mm)",
                                                           "Severely wasted (< -115mm)"))]

############################################################
#' Section 9. Ensure all categorical variables are factors
############################################################
participant_data[, house_fuel_firewood := factor(house_fuel_firewood, c(FALSE, TRUE), c("no", "yes"))]
participant_data[, house_fuel_charcoal := factor(house_fuel_charcoal, c(FALSE, TRUE), c("no", "yes"))]
participant_data[, symptom_respiratory := factor(symptom_respiratory, c(FALSE, TRUE), c("no", "yes"))]
participant_data[, symptom_cough := factor(participant_data$symptom_cough, c(0, 1), c("no", "yes"))]
participant_data[, symptom_sorethroat := factor(participant_data$symptom_sorethroat, c(0, 1), c("no", "yes"))]
participant_data[, symptom_headache := factor(participant_data$symptom_headache, c(0, 1), c("no", "yes"))]
participant_data[, symptom_fever := factor(participant_data$symptom_fever, c(0, 1), c("no", "yes"))]
participant_data[, symptom_diarrhoea := factor(participant_data$symptom_diarrhoea, c(0, 1), c("no", "yes"))]
participant_data[, pneumonia_6m := factor(participant_data$pneumonia_6m, c(0, 1), c("no", "yes"))]
participant_data[, morbidity_sicklecell := factor(participant_data$morbidity_sicklecell, c(0, 1), c("no", "yes"))]
participant_data[, morbidity_asthma := factor(participant_data$morbidity_asthma, c(0, 1), c("no", "yes"))]
participant_data[, morbidity_diabetes := factor(participant_data$morbidity_diabetes, c(0, 1), c("no", "yes"))]
participant_data[, substance_tobacco := factor(participant_data$substance_tobacco, c(0, 1), c("no", "yes"))]
participant_data[, substance_khat := factor(participant_data$substance_khat, c(0, 1), c("no", "yes"))]
participant_data[, antibiotics := factor(participant_data$abs, c(0, 1), c("no", "yes"))]
participant_data[, house_substance_use_smoke := factor(house_substance_use_smoke, c(FALSE, TRUE), c("no", "yes"))]
participant_data[, house_substance_use_snuff := factor(house_substance_use_snuff, c(FALSE, TRUE), c("no", "yes"))]
participant_data[, house_substance_use_khat := factor(house_substance_use_khat, c(FALSE, TRUE), c("no", "yes"))]