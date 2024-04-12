#############################################################################################################################
#' Section 1. Pneumococcal carriage regression
#############################################################################################################################

#' Logistic regression
#' - This analysis assesses the relationship between the collected covariates and NP carriage prevalence
#' - We do not use weighted estimates for the regression, i.e. no post-stratification or finite population corrections
#' - We a priori adjusted estimates for age and sex

#' Check if need to add age as continuous or categorical variable
#' There was strong evidence in favour of more complex model (age as categorical variable)
carriage_age_cont = glm(pneu_carr_final ~ participant_age_y, participant_data[!is.na(participant_age_y) & !is.na(agegrp)], family=stats::binomial(link="logit"))
carriage_age_cat = glm(pneu_carr_final ~ agegrp, participant_data[!is.na(participant_age_y) & !is.na(agegrp)], family=stats::binomial(link="logit"))
lmtest::lrtest(carriage_age_cont, carriage_age_cat)

#' Outcome variable used in the regression
outcome = "pneu_carr_final"

#' A-priori defined confounding variables to add to all regressions in all ages
priori_conf_all = c("agegrp", "participant_sex")

#' Variables to be assessed for all ages
covariates_all = c(
  #household demographics
  "Household size" = "household_size",
  "Household members <5y" = "household_size_u5",
  #shelter
  "House leakage" = "house_leakage",
  "House draft" = "house_draft",
  "Fuel firewood" = "house_fuel_firewood",
  "Fuel charcoal" = "house_fuel_charcoal",
  "Ventilation" = "house_ventilation",
  #health
  "Antibiotic use" = "antibiotics",
  "Respiratory symptoms" = "symptom_respiratory",
  "Cough" = "symptom_cough",
  "Sore throat" = "symptom_sorethroat",
  "Headache" = "symptom_headache",
  "Fever" = "symptom_fever",
  "Diarrhoea" = "symptom_diarrhoea",
  "Pneumonia 6m" = "pneumonia_6m",
  "Sickle Cell" = "morbidity_sicklecell",
  "Asthma" = "morbidity_asthma",
  "Diabetes" = "morbidity_diabetes",
  #substances
  "Tobacco" = "substance_tobacco",
  "Khat" = "substance_khat",
  "Household smoke" = "house_substance_use_smoke",
  "Household snuff" = "house_substance_use_snuff",
  "Household khat" = "house_substance_use_khat",
  #contacts
  "Total number of direct contacts" = "contacts_total_reported",
  "Total number of physical contacts" = "contacts_physical_total_reported")

#' A-priori defined confounding variables to add to all regressions in the <5 subgroup
#' - Note that age in months is included as a continuous variable in the <5 subgroup analyses
priori_conf_u5 = c("participant_age_month", "participant_sex")

#' Variables to be assessed in participants <5 only (anthropometrics were only collected for <5)
covariates_u5 = c(
  "Weight-for-age z-score" = "wfaz",
  "Weight-for-height z-score" = "wfhz",
  "Height-for-age z-score" = "hfaz",
  "MUAC (in cm)" = "muac")

#' Categories used to group models in table
covariates_categories = c(
  "Demographic characteristics" = "Household size",
  "Demographic characteristics" = "Household members <5y",
  "Shelter quality" = "House leakage",
  "Shelter quality" = "House draft",
  "Indoor air pollution" = "Fuel firewood",
  "Indoor air pollution" = "Fuel charcoal",
  "Indoor air pollution" = "Ventilation",
  "Current health" = "Antibiotic use",
  "Current health" = "Respiratory symptoms",
  "Current health" = "Cough",
  "Current health" = "Sore throat",
  "Current health" = "Headache",
  "Current health" = "Fever",
  "Current health" = "Diarrhoea",
  "Morbidities" = "Pneumonia 6m",
  "Morbidities" = "Sickle Cell",
  "Morbidities" = "Asthma",
  "Morbidities" = "Diabetes",
  "Individual substance use" = "Tobacco",
  "Individual substance use" = "Khat",
  "Household substance use" = "Household smoke",
  "Household substance use" = "Household snuff",
  "Household substance use" = "Household khat",
  "Contact behaviour" = "Total number of direct contacts",
  "Contact behaviour" = "Total number of physical contacts",
  "Malnutrition in U5" = "Weight-for-age z-score",
  "Malnutrition in U5" = "Weight-for-height z-score",
  "Malnutrition in U5" = "Height-for-age z-score",
  "Malnutrition in U5" = "MUAC-for-age z-score",
  "Malnutrition in U5" = "MUAC (in cm)")

#' This runs one single logistic regression model for each covariate with the selected outcome (pneumonia)
#' - a priori defined variables are automatically added to each model

#' Analysis in all age-groups
regression_all = 1:length(covariates_all) %>%
  lapply(function(l, covariates_all){
    covars = covariates_all[[l]]
    x = as.formula(paste0(outcome, "~", paste0(c(paste0(priori_conf_all, collapse="+"),
                                                 paste0(covars, collapse="+")), collapse="+"))) %>%
      glm(participant_data, family=stats::binomial(link="logit")) %>% parseGLM(covars) %>% .[, model := l]
    #use LRT for p-value of categorical variables with >2 levels
    if(x[!is.na(pval), .N] > 1){
      x[, pval := NA_real_]
      base_model = as.formula(paste0(outcome, "~", paste0(c(paste0(priori_conf_all, collapse="+")), collapse="+"))) %>%
        glm(participant_data %>% subset(!is.na(get(covars))), family=stats::binomial(link="logit"))
      covar_model = as.formula(paste0(outcome, "~", paste0(c(paste0(priori_conf_all, collapse="+"),
                                                             paste0(covars, collapse="+")), collapse="+"))) %>%
        glm(participant_data, family=stats::binomial(link="logit"))
      lrtp = lmtest::lrtest(base_model, covar_model)
      #default digits_pval = 3 in parseGLM
      x[.N, pval := round(na.omit(lrtp[["Pr(>Chisq)"]])[1], 3)]
    }
    return(x)},
    covariates_all) %>% rbindlist

#' Analysis in the U5
regression_u5 = 1:length(covariates_u5) %>%
  lapply(function(l, covariates_u5){
    covars = covariates_u5[[l]]
    x = as.formula(paste0(outcome, "~", paste0(c(paste0(priori_conf_u5, collapse="+"),
                                                 paste0(covars, collapse="+")), collapse="+"))) %>%
    #x = as.formula(paste0(outcome, "~", paste0(covars, collapse="+"))) %>%
      glm(participant_data_u5, family=stats::binomial(link="logit")) %>%
      parseGLM(covars) %>% .[, model := l]
    return(x)}, covariates_u5) %>% rbindlist

#############################################################################################################################
#' Section 2. Generate Table 2. Association between risk factors and pneumococcal carriage.
#############################################################################################################################
#' Format estimates
regression_all[, variable := factor(variable, covariates_all, names(covariates_all))]
regression_u5[, variable := factor(variable, covariates_u5, names(covariates_u5))]

#' Combine estimates in one table
regression = rbind(regression_all[, data := "all"], regression_u5[, data := "<u5 y"])
regression[, category := factor(variable, covariates_categories, names(covariates_categories))]
regression[, ci_format := paste0(ci_low," - ", ci_high)]
regression[is.na(ci_low) | is.na(ci_high), ci_format := NA_character_]
regression[type == "continuous", c("option", "variable") := .(variable, NA_character_)]
table2_regression_carriage = regression[, c("category", "variable", "option", "est", "ci_format", "pval", "N")]

if(make_table) table2_regression_carriage %>%
  kblOut(booktabs=T, align=c("l", "l","l","l", "l","l", "l"), linesep = "",
         col.names=c("", "Variable", "", "OR", "95% CI", "p-value", "N"), other_functions = list(
           function(x) collapse_rows(x, 1:2, row_group_label_position = "stack",
                                     row_group_label_fonts = list(list(bold=T, italic=T), list(bold=F, italic=T))),
           function(x) kable_styling(x, latex_options = "scale_down")), out_name = "table2_regression_carriage")

#############################################################################################################################
#' Section 3. Pneumococcal density regression
#############################################################################################################################
#' Linear regression
#' - This analysis assesses the relationship between the collected covariates and NP density of all pneumococcal serotypes
#' - We do not use weighted estimates for the regression
#' - We a priori adjusted estimates for age and sex

#' Check if need to add age as continuous or categorical variable
#' There was no evidence in favour of more complex model (age as categorical variable)
#' - however, we kept age as a categorical variable for consistency with our main analysis
density_age_cont = glm(log10(pneu_dens_final) ~ participant_age_y, participant_data[!is.na(participant_age_y) & !is.na(agegrp) & pneu_dens_final > 0], family=stats::gaussian)
density_age_cat = glm(log10(pneu_dens_final) ~ agegrp, participant_data[!is.na(participant_age_y) & !is.na(agegrp) & pneu_dens_final > 0], family=stats::gaussian)
lmtest::lrtest(density_age_cont, density_age_cat)

#' calculate log10 density
participant_data[, lpneudens := log10(pneu_dens_final)]
participant_data_u5[, lpneudens := log10(pneu_dens_final)]

#' Outcome variable to be used in the regression
outcome = "lpneudens"

#' Same as in main regression:
#' priori_conf_all
#' priori_conf_u5
#' covariates_all
#' covariates_u5
#' covariates_categories

#' This runs one single linear regression model for each covariate with the selected outcome
#' - a priori defined variables are automatically added to each model
#' Analysis in all age-groups
regression_all = 1:length(covariates_all) %>%
  lapply(function(l, covariates_all){
    covars = covariates_all[[l]]
    x = as.formula(paste0(outcome, "~", paste0(c(paste0(priori_conf_all, collapse="+"),
                                                 paste0(covars, collapse="+")), collapse="+"))) %>%
      glm(participant_data[pneu_dens_final > 0], family=stats::gaussian) %>% parseGLM(covars, exp = FALSE) %>% .[, model := l]
    #use LRT for p-value of categorical variables with >2 levels
    if(x[!is.na(pval), .N] > 1){
      x[, pval := NA_real_]
      base_model = as.formula(paste0(outcome, "~", paste0(c(paste0(priori_conf_all, collapse="+")), collapse="+"))) %>%
        glm(participant_data[pneu_dens_final > 0] %>% subset(!is.na(get(covars))), family=stats::gaussian)
      covar_model = as.formula(paste0(outcome, "~", paste0(c(paste0(priori_conf_all, collapse="+"),
                                                             paste0(covars, collapse="+")), collapse="+"))) %>%
        glm(participant_data[pneu_dens_final > 0] %>% subset(!is.na(get(covars))), family=stats::gaussian)
      lrtp = lmtest::lrtest(base_model, covar_model)
      #default digits_pval = 3 in parseGLM
      x[.N, pval := round(na.omit(lrtp[["Pr(>Chisq)"]])[1], 3)]
    }
    return(x)},
    covariates_all) %>% rbindlist

#' Analysis in the U5
regression_u5 = 1:length(covariates_u5) %>%
  lapply(function(l, covariates_u5){
    covars = covariates_u5[[l]]
    x = as.formula(paste0(outcome, "~", paste0(c(paste0(priori_conf_u5, collapse="+"),
                                                 paste0(covars, collapse="+")), collapse="+"))) %>%
      #x = as.formula(paste0(outcome, "~", paste0(covars, collapse="+"))) %>%
      glm(participant_data_u5[pneu_dens_final > 0], family=stats::gaussian) %>%
      parseGLM(covars, exp = FALSE) %>% .[, model := l]
    return(x)}, covariates_u5) %>% rbindlist

#############################################################################################################################
#' Section 2. Generate Supplemental Table C3. Association between risk factors and pneumococcal density.
#############################################################################################################################
#' Format estimates
regression_all[, variable := factor(variable, covariates_all, names(covariates_all))]
regression_u5[, variable := factor(variable, covariates_u5, names(covariates_u5))]

#' Combine estimates in one table
regression = rbind(regression_all[, data := "all"], regression_u5[, data := "<u5 y"])
regression[, category := factor(variable, covariates_categories, names(covariates_categories))]
regression[, ci_format := paste0(ci_low," - ", ci_high)]
regression[is.na(ci_low) | is.na(ci_high), ci_format := NA_character_]
regression[type == "continuous", c("option", "variable") := .(variable, NA_character_)]
tableSC3_regression_density = regression[, c("category", "variable", "option", "est", "ci_format", "pval", "N")]

if(make_table) tableSC3_regression_density %>%
  kblOut(booktabs=T, align=c("l", "l","l","l", "l","l", "l"), linesep = "",
         col.names=c("", "Variable", "", "Mean log-difference", "95% CI", "p-value", "N"), other_functions = list(
           function(x) collapse_rows(x, 1:2, row_group_label_position = "stack",
                                     row_group_label_fonts = list(list(bold=T, italic=T), list(bold=F, italic=T))),
           function(x) kable_styling(x, latex_options = "scale_down")), out_name = "tableSC3_regression_density")

