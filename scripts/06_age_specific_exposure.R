#############################################################################################################################
#' Section 1. Get Digaale contact data from Zenodo
#' - Nb, this is the same as the data in the kevinvzandvoort/espicc-somaliland-digaale-survey-2019 repository
#############################################################################################################################
#' Get data for socialmixr
digaale_contact_data = socialmixr::get_survey("https://zenodo.org/record/5226280")

#' The estimated population size in Digaale (for provided age groups)
digaale_survey_population =
  data.table::fread("https://zenodo.org/record/7071876/files/espicc_somaliland_digaale_survey_population.csv")

#' Note that weekends fall on Fridays and Saturdays in Somaliland.
#' - The dayofweek variable provided in the dataset has been kept
#'   consistent with R defaults (0: Sunday to 6: Saturday)
digaale_contact_data$participants[, c("dayofweek", "dayofweek_name", "weekend")] %>%
  unique %>% setorder(dayofweek) %>% .[]
digaale_contact_data$participants[, dayofweek := ifelse(dayofweek == 6, 0, dayofweek + 1)]

#' The contact matrix can then be constructed as follows
#' - The provided survey_population can be used to construct a population representative matrix for Digaale IDP camp
#' - As the sample is not self-weighing (oversampling of young age groups), we apply the survey_weight as weights
digaale_survey_population = poststratification_strata_list$age
digaale_survey_population[, lower.age.limit := c(0, 2, 6, 15, 30, 50)]
digaale_survey_population = digaale_survey_population[, c("lower.age.limit", "Freq")]
colnames(digaale_survey_population) = c("lower.age.limit", "population")

#' calculate contact matrix, adjusted for reciprocity of contacts
digaale_contact_matrix = digaale_contact_data %>%
  socialmixr::contact_matrix(survey.pop = digaale_survey_population,
                             age.limits = digaale_survey_population$lower.age.limit,
                             symmetric = TRUE, weights = "survey_weight", weigh.dayofweek = TRUE)

#############################################################################################################################
#' Section 2. Combine with prevalence by age estimates
#############################################################################################################################
prev_by_age = svyMean2(~pneu_carr_final, participant_data_design_ps, na.rm=T, by="agegrp") %>%
  setorder(agegrp) %>% .[, mean]
carrier_contacts = t(digaale_contact_matrix$matrix) %*% diag(prev_by_age)
rel_carrier_contacts = carrier_contacts/rowSums(carrier_contacts)
rel_carrier_contacts = as.data.table(rel_carrier_contacts)
colnames(rel_carrier_contacts) = age_groups_sample$name
rel_carrier_contacts[, contactor_age := age_groups_sample$name]

#' Note that socialmixr contactors are in rows, contactees are in columns
cm = t(digaale_contact_matrix$matrix)
carrier_contacts = t(t(cm) %*% diag(prev_by_age))

#' scale so columns sum to 1
rel_carrier_contacts = carrier_contacts %*% diag(1/colSums(carrier_contacts))
rel_carrier_contacts = as.data.table(rel_carrier_contacts)
colnames(rel_carrier_contacts) = age_groups_sample$name
rel_carrier_contacts[, contactee_age := age_groups_sample$name]

#############################################################################################################################
#' Section 3. Generate Figure 4: The contribution of different age groups towards the age specific exposure to pneumococcus.
#############################################################################################################################
figure4_proportion_infected_contacts = rel_carrier_contacts %>%
  melt(id.vars = "contactee_age", variable.name = "contactor_age", value.name = "rel_transmission") %>%
  .[, contactor_age := factor(contactor_age, age_groups_sample$name)] %>%
  .[, contactee_age := factor(contactee_age, age_groups_sample$name)] %>%
  ggplot(aes(x=contactor_age, y=rel_transmission, fill=factor(contactee_age, rev(age_groups_sample$name))))+
  geom_col()+
  scale_fill_brewer(palette="Greens")+
  scale_y_continuous(labels = scales::percent)+
  labs(x="Contactor age group", fill="Contactee age group", y="Proportion of infected contacts")+
  theme_minimal()+
  theme_multiplot+
  theme(legend.position = "right", axis.text.x = element_text(angle=45, hjust=1))

for(ext in c("png", "pdf", "tiff", "eps"))
  ggsave(sprintf("%s/output/%s/figures/%s/figure4_proportion_age_exposure_carriers.%s", analysis_dir, OUTPUT_DIR, ext, ext),
         plot = figure4_proportion_infected_contacts, width = plot_single_col*1.25, height = plot_single_col*0.75, units = "in", dpi = 300)

#############################################################################################################################
#' Section 4. Bootstrap analysis to assess uncertainty around estimates.
#############################################################################################################################
set.seed(123)
boot_samples = 1000

#' bootstrap contact matrices
digaale_contact_boot_matrices =
  suppressMessages(replicate(n = boot_samples,
            socialmixr::contact_matrix(digaale_contact_data, survey.pop = digaale_survey_population,
                                       age.limits = digaale_survey_population$lower.age.limit,
                                       symmetric = TRUE, weights = "survey_weight", weigh.dayofweek = TRUE,
                                       sample.participants = TRUE)))
digaale_contact_boot_matrices = digaale_contact_boot_matrices[seq(1, 3*boot_samples, by=3)]
  
#' Need to take new bootstrap samples for carriage estimates
boot_data_init = participant_data[participant_consent == "yes" & !is.na(fpc)]
digaale_prevalence_boot = c(1:boot_samples) %>% lapply(function(i){
  #bootstrapped sample
  boot_pop = boot_data_init[sample(1:.N, .N, TRUE)]
  
  #make Survey object
  boot_participant_data_design_svy =
    survey::svydesign(ids=~1, probs=NULL, strata=~stype, fpc=~fpc, weights=~pw,
                      data=boot_pop)  
  
  #poststratify
  boot_participant_data_design_ps = boot_participant_data_design_svy %>%
    postStratify(~stype, poststratification_strata, partial=TRUE)
  
  #trim weights
  wlow = quantile(weights(boot_participant_data_design_ps), 0.05)
  whigh = quantile(weights(boot_participant_data_design_ps), 0.95)
  boot_participant_data_design_ps = boot_participant_data_design_ps %>%
    survey::trimWeights(whigh, wlow)
  
  #get bootstrapped prevalence by age
  prev_by_age = svyMean2(~pneu_carr_final, boot_participant_data_design_ps, na.rm=T, by="agegrp") %>%
    setorder(agegrp) %>%
    .[, mean]
  
  return(prev_by_age)
})

#' Combine contact and carriage estimates in each bootstrapped sample
bootstrapped_estimates = c(1:boot_samples) %>% lapply(function(i) {
  cm = t(digaale_contact_boot_matrices[[i]])
  carrier_contacts = t(t(cm) %*% diag(digaale_prevalence_boot[[i]]))
  rel_carrier_contacts = carrier_contacts %*% diag(1/colSums(carrier_contacts))
  rel_carrier_contacts = as.data.table(rel_carrier_contacts)
  colnames(rel_carrier_contacts) = age_groups_sample$name
  rel_carrier_contacts[, contactee_age := age_groups_sample$name]
  
  rel_carrier_contacts =
    rel_carrier_contacts %>% melt(id.vars="contactee_age", variable.name = "contactor_age") %>%
    .[, b := i]
  
  return(rel_carrier_contacts)
})
bootstrapped_estimates = bootstrapped_estimates %>% rbindlist

#' calculate summary estimates
bootstrapped_estimates_summary = bootstrapped_estimates[, .(mean = mean(value), median = median(value),
                           low95 = quantile(value, 0.025),
                           high95 = quantile(value, 0.975),
                           low50 = quantile(value, 0.25),
                           high50 = quantile(value, 0.75)), by=c("contactor_age", "contactee_age")]

#' format for table
bootstrapped_estimates_summary[, val := sprintf("%s%% (%s - %s)", round(median*100, 1), round(low95*100, 1), round(high95*100, 1))]
bootstrapped_estimates_summary = bootstrapped_estimates_summary %>% dcast(contactee_age~contactor_age, value.var = "val") %>%
  .[, contactee_age := factor(contactee_age, age_groups_sample$name)] %>%
  setorder(contactee_age) %>% .[]
bootstrapped_estimates_summary[, c("contactee_age", age_groups_sample$name), with=FALSE]

#############################################################################################################################
#' Section 4. Generate Supplemental Table C6. The contribution of different age groups towards the age specific exposure
#'            to pneumococcus.
#############################################################################################################################
if(make_table) bootstrapped_estimates_summary %>%
  kblOut(booktabs=T, align=c("l", "l","l","l", "l","l", "l"), linesep = "",
         col.names=c("Contactee age", "[0, 2)", "[2, 6)", "[6, 15)", "[15, 29)", "[30, 49)", "50+"),
         other_functions = list(
           function(x) kable_styling(x, latex_options = "scale_down")), out_name = "tableSC6_bootstrapped_WAIFW")
