#' Some values in this script are read from other files. These values are used to calculate the post-stratification weights.
#' They are equivalent to those in the https://github.com/kevinvzandvoort/espicc-somaliland-digaale-survey-2019 repository.
#' Download the kevinvzandvoort/espicc-somaliland-digaale-survey-2019 repository to see calculations for these values.

#' household grouping used when poststratificatied by household size
household_size_groups = readRDS("./data/poststratification_data/household_size_groups.RDS")
household_size_groups_u5 = readRDS("./data/poststratification_data/household_size_groups_u5.RDS")

#' Post-stratification is done on the joint age- and sex-distribution of the population
#'  - The age- and sex-distribution is calculated as the observed distribution in all households included in the survey
#'    inflated to account for households not included in the sample
#'  - We also conduct several sensitivity analyses based on the weights that will be used
poststratification_strata_list = readRDS("./data/poststratification_data/poststratification_strata_list.RDS")
poststratification_strata_u5_list = readRDS("./data/poststratification_data/poststratification_strata_u5_list.RDS")

#' Retrieve the population distribution for the poststratification strata that are used
poststratification_strata = poststratification_strata_list[[POSTSTRATIFICATION_STRATA]]
poststratification_strata_u5 = poststratification_strata_u5_list[[POSTSTRATIFICATION_STRATA]] %>% .[, c("stype", "Freq")]

############################################################
#' Section 1. Participant data
############################################################  
#' We apply the post-stratification strata to our participant using the appropriate methods from the Survey package
#' - the population size in each age- and sex-stratum is used to apply a finite population correction (fpc)
participant_data = participant_data %>%
  .[, pw := 1] %>%
  .[, pw := pw * .N/sum(pw)] %>%
  .[, household_size_group := factor(as.character(household_size), household_size_groups$values,
                                     household_size_groups$names)] %>%
  .[, stype_age_sex_householdsize := paste(household_size_group, participant_age_group_sample, participant_sex, sep="-")] %>%
  .[, stype_age_householdsize := paste(household_size_group, participant_age_group_sample, sep="-")] %>%
  .[, stype_age_sex := paste(participant_sex, participant_age_group_sample, sep="-")] %>%
  .[, stype_age := participant_age_group_sample] %>%
  .[, stype_none := "none"]

participant_data = participant_data %>%
  .[, stype := get(paste0("stype_", POSTSTRATIFICATION_STRATA))] %>%
  merge(poststratification_strata, by="stype", all.x = TRUE) %>% .[, fpc := Freq]

#' Create a survey design object to process participant estimates
#' - we set partial=TRUE as there are no male participants aged 30-49yo who live in household of size 5-6 or 7-12
#'   in our sample
participant_data_design_svy = survey::svydesign(ids=~1, probs=NULL, strata=~stype, fpc=~fpc, weights=~pw,
                                               data=participant_data[participant_consent == "yes" & !is.na(fpc)])
participant_data_design_ps = participant_data_design_svy %>% postStratify(~stype, poststratification_strata, partial=TRUE)

#' Trim weights that are outliers so weights don't overinfluence data
wlow = quantile(weights(participant_data_design_ps), 0.05)
whigh = quantile(weights(participant_data_design_ps), 0.95)
participant_data_design_ps = participant_data_design_ps %>%
  survey::trimWeights(whigh, wlow)

####################################################################
#' Section 2. Participant data <5 (for nutrition data subanalyses)
####################################################################
participant_data_u5 = 
  participant_data[participant_age_month < 5*12, -c("stype_age_sex", "stype_age", "stype_none", "Freq", "fpc")]

participant_data_u5 %>%
  .[, household_size_group_u5 := factor(as.character(household_size),
                                        household_size_groups_u5$values, household_size_groups_u5$names)] %>%
  .[, participant_age_years := floor(participant_age_month/12)] %>%
  .[, stype_age_sex_householdsize := paste(household_size_group_u5, participant_age_years, participant_sex, sep="-")] %>%
  .[, stype_age_householdsize := paste(household_size_group_u5, participant_age_years, sep="-")] %>%
  .[, stype_age_sex := paste(participant_sex, participant_age_years, sep="-")] %>%
  .[, stype_age := participant_age_years] %>%
  .[, stype_none := "none"]

#' Set correct stype
participant_data_u5 = participant_data_u5 %>%
  .[, stype := get(paste0("stype_", POSTSTRATIFICATION_STRATA))] %>%
  merge(poststratification_strata_u5, by="stype") %>% .[, fpc := Freq] %>% .[, -"Freq"]

participant_data_u5_design_svy = survey::svydesign(ids=~1, probs=NULL, strata=~stype, fpc=~fpc, weights=~pw,
                                                  data=participant_data_u5[participant_consent == "yes" & !is.na(fpc)])

participant_data_u5_design_ps = participant_data_u5_design_svy %>%
  postStratify(~stype, poststratification_strata_u5, partial=TRUE)

#' Trim weights to remove any outliers
wlow = quantile(weights(participant_data_u5_design_ps), 0.05)
whigh = quantile(weights(participant_data_u5_design_ps), 0.95)
participant_data_u5_design_ps = participant_data_u5_design_ps %>% trimWeights(whigh, wlow)

####################################################################
#' Section 3. Pilot data
#' - used to split datasets in those with/without pilot data
####################################################################
participant_data = participant_data %>% merge(lab_data[, c("pid", "pilot")], by="pid", all.x = TRUE)
participant_data_u5 = participant_data_u5 %>% merge(lab_data[, c("pid", "pilot")], by="pid", all.x = TRUE)

participant_data_pilot_design_svy =
  survey::svydesign(ids=~1, probs=NULL, strata=~stype, fpc=~fpc, weights=~pw,
                    data=participant_data[participant_consent == "yes" & !is.na(fpc) & pilot == TRUE])
participant_data_pilot_design_ps = participant_data_pilot_design_svy %>%
  postStratify(~stype, poststratification_strata, partial=TRUE)

#' Trim weights that are outliers so weights don't overinfluence data
wlow = quantile(weights(participant_data_pilot_design_ps), 0.05)
whigh = quantile(weights(participant_data_pilot_design_ps), 0.95)
participant_data_pilot_design_ps = participant_data_pilot_design_ps %>%
  survey::trimWeights(whigh, wlow)

participant_data_u5_pilot_design_svy =
  survey::svydesign(ids=~1, probs=NULL, strata=~stype, fpc=~fpc, weights=~pw,
                    data=participant_data_u5[participant_consent == "yes" & !is.na(fpc) & pilot == TRUE])

participant_data_u5_pilot_design_ps = participant_data_u5_pilot_design_svy %>%
  postStratify(~stype, poststratification_strata_u5, partial=TRUE)

#' Trim weights to remove any outliers
wlow = quantile(weights(participant_data_u5_pilot_design_ps), 0.05)
whigh = quantile(weights(participant_data_u5_pilot_design_ps), 0.95)
participant_data_u5_pilot_design_ps = participant_data_u5_pilot_design_ps %>% trimWeights(whigh, wlow)

#' make alternative datasets without pilot data
participant_data_nopilot_design_svy =
  survey::svydesign(ids=~1, probs=NULL, strata=~stype, fpc=~fpc, weights=~pw,
                    data=participant_data[participant_consent == "yes" & !is.na(fpc) & pilot == FALSE])
participant_data_nopilot_design_ps = participant_data_nopilot_design_svy %>%
  postStratify(~stype, poststratification_strata, partial=TRUE)

#' Trim weights that are outliers so weights don't overinfluence data
wlow = quantile(weights(participant_data_nopilot_design_ps), 0.05)
whigh = quantile(weights(participant_data_nopilot_design_ps), 0.95)
participant_data_nopilot_design_ps = participant_data_nopilot_design_ps %>%
  survey::trimWeights(whigh, wlow)

participant_data_u5_nopilot_design_svy =
  survey::svydesign(ids=~1, probs=NULL, strata=~stype, fpc=~fpc, weights=~pw,
                    data=participant_data_u5[participant_consent == "yes" & !is.na(fpc) & pilot == FALSE])

participant_data_u5_nopilot_design_ps = participant_data_u5_nopilot_design_svy %>%
  postStratify(~stype, poststratification_strata_u5, partial=TRUE)

#' Trim weights to remove any outliers
wlow = quantile(weights(participant_data_u5_nopilot_design_ps), 0.05)
whigh = quantile(weights(participant_data_u5_nopilot_design_ps), 0.95)
participant_data_u5_nopilot_design_ps = participant_data_u5_nopilot_design_ps %>% trimWeights(whigh, wlow)