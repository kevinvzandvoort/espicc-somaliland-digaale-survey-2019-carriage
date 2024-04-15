#' Data from https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009389
#' - The median estimates, and lower and upper bounds of the 95% credibility interval, are listed for each serotype.
#' - Invasiveness is defined as the progression rate from carriage to invasive disease, in cases per carrier per year

#############################################################################################################################
#' Section 1. Compare data serotype match with Digaale
#############################################################################################################################
#' check data match for children
circulating_children = lab_data %>% merge(participant_data[, c("participant_age_y","pid")], all.x=TRUE) %>%
  .[participant_age_y < 18 & pneu_carr_final == 1,
    c(sero1_final, sero2_final, sero3_final, sero4_final, sero5_final)] %>% unlist %>% unique()
circulating_children = na.omit(circulating_children) %>% subset(. != "")
data_children = readxl::read_excel("./data/other_data/pcbi.1009389.s048_children.xlsx")
data_children = as.data.table(data_children)
length(circulating_children %>% subset(. %in% data_children$Serotype))
length(circulating_children %>% subset(!. %in% data_children$Serotype))
circulating_children %>% subset(!. %in% data_children$Serotype)
#' 20B not present - 20 is present, use 20 (as it is the same serogroup)
#' NT3b, NT4a, NT2, and NT4b not present, use average over all values OR assume invasiveness is 0

#' check data match for adults
circulating_adults = lab_data %>% merge(participant_data[, c("participant_age_y","pid")], all.x=TRUE) %>%
  .[participant_age_y >= 18 & pneu_carr_final == 1,
    c(sero1_final, sero2_final, sero3_final, sero4_final, sero5_final)] %>% unlist %>% unique()
circulating_adults = na.omit(circulating_adults) %>% subset(. != "")
data_adults = readxl::read_excel("./data/other_data/pcbi.1009389.s049_adults.xlsx")
data_adults = as.data.table(data_adults)
length(circulating_adults %>% subset(. %in% data$Serotype))
circulating_adults %>% subset(!. %in% data$Serotype)
#' 20B not present - 20 is present, use 20 (as it is the same serogroup)
#' 41A is not present, average over all values
#' NT3b, NT2, and NT4a are not present, average over all values OR assume invasiveness is 0 

#############################################################################################################################
#' Section 2. Refit Lochen et al estimates to get distributions
#############################################################################################################################
optDist = function(x, i=1){
  log(sum((data_children[i, c(invasiveness_lower, invasiveness_upper)] - qlnorm(c(0.025, 0.975), data_children[i, log(invasiveness)], x))^2)  )
}
for(i in 1:nrow(data_children)){
  fit = optimize(optDist, c(1e-12, 1e2), c(i=i))
  data_children[i, c("ln_mu", "ln_sd", "fit_obj") := .(log(invasiveness), fit$minimum, fit$objective)]
}
invasiveness_children = data_children

optDist = function(x, i=1){
  log(sum((data_adults[i, c(invasiveness_lower, invasiveness_upper)] - qlnorm(c(0.025, 0.975), data_adults[i, log(invasiveness)], x))^2)  )
}
for(i in 1:nrow(data_adults)){
  fit = optimize(optDist, c(1e-12, 1e2), c(i=i))
  data_adults[i, c("ln_mu", "ln_sd", "fit_obj") := .(log(invasiveness), fit$minimum, fit$objective)]
}
invasiveness_adults = data_adults

#############################################################################################################################
#' Section 3. Generate Supplemental Figure D1: Age and serotype specific invasiveness values.
#############################################################################################################################
data_plot = rbind(data_children[, age := "children"], data_adults[, age := "adults"])
data_plot[, age := factor(age, c("children", "adults"))]
setorder(data_plot, age, fit_obj)

figure_ipd_refitted = data_plot %>% copy %>% .[, type:="posterior"] %>%
  rbind(copy(.) %>% .[, c("invasiveness", "invasiveness_lower", "invasiveness_upper", "type") :=
                        .(qlnorm(0.5, ln_mu, ln_sd), qlnorm(0.025, ln_mu, ln_sd), qlnorm(0.975, ln_mu, ln_sd), "fitted")]) %>%
  ggplot(aes(x=factor(paste0(age, Serotype),
                      data_plot[, paste0(age, Serotype)],
                      data_plot[, Serotype]),
             y=invasiveness, ymin=invasiveness_lower, ymax=invasiveness_upper,
             colour=factor(type, c("fitted", "posterior"))))+
  facet_grid(.~age, scales = "free")+
  geom_errorbar(position=position_dodge(width=1))+
  geom_point(position=position_dodge(width=1))+
  scale_y_log10()+
  scale_colour_manual(values=c("posterior" = "#0D5257", "fitted" = "#00BF6F"))+
  theme_bw()+
  labs(x="Serotype", colour="Distribution")+
  coord_flip()+
  theme(legend.position="bottom")

#ggsave("./refitted_distributions.png", width=6.954, height=10*(6.954/8))
for(ext in c("png", "pdf", "tiff", "eps"))
  ggsave(sprintf("%s/output/%s/figures/%s/figureSD1_refitted_invasiveness_estimates.%s", analysis_dir, OUTPUT_DIR, ext, ext),
         plot = figure_ipd_refitted, width = (page_width - page_outer_margin*2)*1.15, height = 10*(6.954/8)*0.9, units = "in", dpi = 300)

#############################################################################################################################
#' Section 4. Take bootstrap samples from datasets
#############################################################################################################################
set.seed(123)
boot_samples = 1000

#' bootstrap population for carriage estimates
boot_data_init = participant_data[participant_consent == "yes" & !is.na(fpc)]
boot_data_init = boot_data_init %>% merge(lab_data, by="pid")

boot_results = digaale_boot_populations = c(1:boot_samples) %>%
  lapply(function(b, boot_samples){
    message(sprintf("%s/%s (%s%%)", b, boot_samples, round(b/boot_samples * 100)))
    
    #' a. create survey design object and poststratify
    boot_pop = boot_data_init[sample(1:.N, .N, TRUE)]
    boot_pop_u5 = copy(boot_pop[participant_age_y < 5])
    boot_pop_u5 %>%
      .[, household_size_group_u5 := factor(as.character(household_size),
                                            household_size_groups_u5$values, household_size_groups_u5$names)] %>%
      .[, participant_age_years := floor(participant_age_month/12)] %>%
      .[, stype_age_sex_householdsize := paste(household_size_group_u5, participant_age_years, participant_sex, sep="-")] %>%
      .[, stype_age_householdsize := paste(household_size_group_u5, participant_age_years, sep="-")] %>%
      .[, stype_age_sex := paste(participant_sex, participant_age_years, sep="-")] %>%
      .[, stype_age := participant_age_years] %>%
      .[, stype_none := "none"]
    boot_pop_u5 = boot_pop_u5[, -"Freq"] %>%
      .[, stype := get(paste0("stype_", POSTSTRATIFICATION_STRATA))]
    
    #' Make sure there is at least one participant in each age group
    while(boot_pop[, .N, by="stype"] %>% .[, any(N < 2)] | boot_pop_u5[, .N, by="stype"] %>% .[, any(N < 2)]){
      boot_pop = boot_data_init[sample(1:.N, .N, TRUE)]
    }
    
    #' b. sample invasiveness value for each circulating serotype, in children and adults
    boot_pop[, invasiveness_novac := 0]
    for(vaccine in c("pneumosil", "pcv10", "pcv13", "pcv15", "pcv20")){
      boot_pop[, sprintf("invasiveness_%s", vaccine) := 0]
    }
    
    invasiveness_children[, val := rlnorm(1, ln_mu, ln_sd), by="Serotype"]
    invasiveness_adults[, val := rlnorm(1, ln_mu, ln_sd), by="Serotype"]
    
    for(s in circulating_children){
      if(s %in% c("NT3b", "NT4a", "NT2", "NT4b")){
        #use 0 - assuming no invasiveness from NEps
        inv = 0
      } else if(s == "20B"){
        inv = invasiveness_children[Serotype == "20", val]
      } else {
        inv = invasiveness_children[Serotype == s, val]
      }
      
      #' no-vaccination invasiveness
      boot_pop[participant_age_y < 18 &
                 (sero1_final == s | sero2_final == s | sero3_final == s | sero4_final == s | sero5_final == s), invasiveness_novac :=
                 invasiveness_novac + inv]
      
      #' invasiveness without serotypes in vaccine
      for(vaccine in c("pneumosil", "pcv10", "pcv13", "pcv15", "pcv20")){
        if(s %in% serotypes_by_vaccine[get(vaccine) == 1, serotype]){
          boot_pop[participant_age_y < 18 &
                     (sero1_final == s | sero2_final == s | sero3_final == s | sero4_final == s | sero5_final == s), sprintf("invasiveness_%s", vaccine) := get(sprintf("invasiveness_%s", vaccine)) + 0]
        } else {
          boot_pop[participant_age_y < 18 &
                     (sero1_final == s | sero2_final == s | sero3_final == s | sero4_final == s | sero5_final == s), sprintf("invasiveness_%s", vaccine) :=
                     get(sprintf("invasiveness_%s", vaccine)) + inv]
        }
      }
    }
    for(s in circulating_adults){
      if(s %in% c("NT3b", "NT4a", "NT2")){
        #use 0 or mean of all values
        #inv = invasiveness_adults[, median(val)]
        inv = 0
      } else if(s == "20B"){
        inv = invasiveness_adults[Serotype == "20", val]
      } else if(s == "19B"){
        inv = invasiveness_adults[Serotype %in% c("19A", "19F"), median(val)]
      } else if(s == "41A"){
        inv = invasiveness_adults[, median(val)]
      } else {
        inv = invasiveness_adults[Serotype == s, val]
      }
      
      #' no-vaccination invasiveness
      boot_pop[participant_age_y >= 18 &
                 (sero1_final == s | sero2_final == s | sero3_final == s | sero4_final == s | sero5_final == s), invasiveness_novac :=
                 invasiveness_novac + inv]
      #' invasiveness without serotypes in vaccine
      for(vaccine in c("pneumosil", "pcv10", "pcv13", "pcv15", "pcv20")){
        if(s %in% serotypes_by_vaccine[get(vaccine) == 1, serotype]){
          boot_pop[participant_age_y >= 18 &
                     (sero1_final == s | sero2_final == s | sero3_final == s | sero4_final == s | sero5_final == s), sprintf("invasiveness_%s", vaccine) := get(sprintf("invasiveness_%s", vaccine)) + 0]
        } else {
          boot_pop[participant_age_y >= 18 &
                     (sero1_final == s | sero2_final == s | sero3_final == s | sero4_final == s | sero5_final == s), sprintf("invasiveness_%s", vaccine) :=
                     get(sprintf("invasiveness_%s", vaccine)) + inv]
        }
      }
    }
    
    #' make ps design for all ages and <5y
    #' design for all age groups
    boot_participant_data_design_svy =
      survey::svydesign(ids=~1, probs=NULL, strata=~stype, fpc=~fpc, weights=~pw,
                        data=boot_pop)  
    
    boot_participant_data_design_ps = boot_participant_data_design_svy %>%
      postStratify(~stype, poststratification_strata, partial=TRUE)
    
    wlow = quantile(weights(boot_participant_data_design_ps), 0.05)
    whigh = quantile(weights(boot_participant_data_design_ps), 0.95)
    boot_participant_data_design_ps = boot_participant_data_design_ps %>%
      survey::trimWeights(whigh, wlow)
    
    #' design for <5y
    boot_pop_u5 = copy(boot_pop[participant_age_y < 5])
    boot_pop_u5 %>%
      .[, household_size_group_u5 := factor(as.character(household_size),
                                            household_size_groups_u5$values, household_size_groups_u5$names)] %>%
      .[, participant_age_years := floor(participant_age_month/12)] %>%
      .[, stype_age_sex_householdsize := paste(household_size_group_u5, participant_age_years, participant_sex, sep="-")] %>%
      .[, stype_age_householdsize := paste(household_size_group_u5, participant_age_years, sep="-")] %>%
      .[, stype_age_sex := paste(participant_sex, participant_age_years, sep="-")] %>%
      .[, stype_age := participant_age_years] %>%
      .[, stype_none := "none"]
    
    boot_pop_u5 = boot_pop_u5[, -"Freq"] %>%
      .[, stype := get(paste0("stype_", POSTSTRATIFICATION_STRATA))] %>%
      merge(poststratification_strata_u5, by="stype") %>% .[, fpc := Freq] %>% .[, -"Freq"]
    
    boot_participant_data_u5_design_svy = survey::svydesign(ids=~1, probs=NULL, strata=~stype, fpc=~fpc, weights=~pw,
                                                            data=boot_pop_u5[participant_consent == "yes" & !is.na(fpc)])
    
    boot_participant_data_u5_design_ps = boot_participant_data_u5_design_svy %>%
      postStratify(~stype, poststratification_strata_u5, partial=TRUE)
    
    #' Trim weights to remove any outliers
    wlow = quantile(weights(boot_participant_data_u5_design_ps), 0.05)
    whigh = quantile(weights(boot_participant_data_u5_design_ps), 0.95)
    boot_participant_data_u5_design_ps = boot_participant_data_u5_design_ps %>% trimWeights(whigh, wlow)
    
    #' calculate difference with vaccination
    estimates = lapply(c("pneumosil", "pcv10", "pcv13", "pcv15", "pcv20"), function(vaccine){
      #sample
      difference_sample = 1 - (boot_participant_data_design_ps$variables[, sum(get(sprintf("invasiveness_%s", vaccine)))]/boot_participant_data_design_ps$variables[, sum(invasiveness_novac)])
      #all ages
      difference_pop_allages = 1 - (svyTotal2(~get(sprintf("invasiveness_%s", vaccine)), boot_participant_data_design_ps, na.rm=T)[, total]/svyTotal2(~invasiveness_novac, boot_participant_data_design_ps, na.rm=T)[, total])
      #under 5y
      difference_pop_u5 = 1 - (svyTotal2(~get(sprintf("invasiveness_%s", vaccine)), boot_participant_data_u5_design_ps, na.rm=T)[, total]/svyTotal2(~invasiveness_novac, boot_participant_data_u5_design_ps, na.rm=T)[, total])
      #by agegroup
      difference_pop_byage = svyTotal2(~get(sprintf("invasiveness_%s", vaccine)), boot_participant_data_design_ps, na.rm=T, by="agegrp") %>%
        setorder(agegrp) %>% .[, c("agegrp", "total")] %>%
        .[, vaccinated := total] %>% .[, -"total"] %>%
        cbind(svyTotal2(~invasiveness_novac, boot_participant_data_design_ps, na.rm=T, by="agegrp") %>%
                setorder(agegrp) %>% .[, total]) %>%
        .[, .(agegrp = agegrp, difference = 1 - vaccinated/V2)]
      return(list(sample = data.table(difference = difference_sample, vaccine = vaccine),
                  pop_all_ages = data.table(difference = difference_pop_allages, vaccine = vaccine),
                  pop_u5 = data.table(difference=difference_pop_u5, vaccine = vaccine),
                  pop_byage = difference_pop_byage[, vaccine := vaccine]))})
    
    return(list(sample = rbindlist(estimates %>% lapply("[[", "sample")) %>% .[, boot_sample := b],
                pop_all_ages = rbindlist(estimates %>% lapply("[[", "pop_all_ages")) %>% .[, boot_sample := b],
                pop_u5 = rbindlist(estimates %>% lapply("[[", "pop_u5")) %>% .[, boot_sample := b],
                pop_byage = rbindlist(estimates %>% lapply("[[", "pop_byage")) %>% .[, boot_sample := b]))
  }, boot_samples)

#############################################################################################################################
#' Section 5. Calculate population-level estimates
#############################################################################################################################
sample_estimates = rbindlist(boot_results %>% lapply("[[", "sample")) %>%
  .[, .(med=(median(difference)), low=(quantile(difference, 0.025)), high=(quantile(difference, 0.975))), by="vaccine"] %>%
  .[, vaccine := factor(vaccine, c("pneumosil", "pcv10", "pcv13", "pcv15", "pcv20"),
                        c("PNEUMOSIL", "Synflorix", "Prevenar 13", "Vaxneuvance", "Prevenar 20"))]

population_estimates_all = rbindlist(boot_results %>% lapply("[[", "pop_all_ages")) %>%
  .[, .(med=median(difference), low=quantile(difference, 0.025), high=quantile(difference, 0.975)), by="vaccine"] %>%
  .[, vaccine := factor(vaccine, c("pneumosil", "pcv10", "pcv13", "pcv15", "pcv20"),
                        c("PNEUMOSIL", "Synflorix", "Prevenar 13", "Vaxneuvance", "Prevenar 20"))]

population_estimates_u5 = rbindlist(boot_results %>% lapply("[[", "pop_u5")) %>%
  .[, .(med=median(difference), low=quantile(difference, 0.025), high=quantile(difference, 0.975)), by="vaccine"] %>%
  .[, vaccine := factor(vaccine, c("pneumosil", "pcv10", "pcv13", "pcv15", "pcv20"),
                        c("PNEUMOSIL", "Synflorix", "Prevenar 13", "Vaxneuvance", "Prevenar 20"))]

population_estimates_byage = rbindlist(boot_results %>% lapply("[[", "pop_byage")) %>%
  .[, .(med=median(difference, na.rm=TRUE), low=quantile(difference, 0.025, na.rm=TRUE), high=quantile(difference, 0.975, na.rm=TRUE)), by=c("vaccine", "agegrp")] %>%
  .[, vaccine := factor(vaccine, c("pneumosil", "pcv10", "pcv13", "pcv15", "pcv20"),
                        c("PNEUMOSIL", "Synflorix", "Prevenar 13", "Vaxneuvance", "Prevenar 20"))]

#############################################################################################################################
#' Section 6. Generate Supplemental Figure D2: Estimated proportion of IPD cases caused by serotypes covered by PCVs.
#############################################################################################################################
figure_ipd_serotypes = population_estimates_byage %>%
  rbind(population_estimates_all %>% .[, agegrp := "All"]) %>%
  .[, agegrp := factor(agegrp, c("All", age_groups_sample$name))] %>%
  setorder(agegrp, vaccine) %>%
  ggplot(aes(x=agegrp, y=med, ymin=low, ymax=high, colour=vaccine))+
  geom_errorbar(width=0.1, position=position_dodge(width=0.5))+
  geom_point(position=position_dodge(width=0.5))+
  scale_y_continuous(limits=c(0, 1), labels=scales::percent)+
  theme_bw()+
  labs(x="Age group", colour="PCV", y="IPD caused by\nserotypes covered by PCV")+
  scale_colour_manual(values=c("PNEUMOSIL" = lshtm_colours$blue,
                               "Synflorix" = lshtm_colours$purple,
                               "Prevenar 13" = lshtm_colours$pink,
                               "Vaxneuvance" = lshtm_colours$lightgreen,
                               "Prevenar 20" = lshtm_colours$yellow))
for(ext in c("png", "pdf", "tiff", "eps"))
  ggsave(sprintf("%s/output/%s/figures/%s/figureSD2_ipd_covered_pcv.%s", analysis_dir, OUTPUT_DIR, ext, ext),
         plot = figure_ipd_serotypes, width = plot_double_col, height = plot_single_col, units = "in", dpi = 300)

#############################################################################################################################
#' Section 7. Generate Table 1b: Proportion of current IPD covered by PCVs
#############################################################################################################################
table1_ipd_pcv = data.table(heading = "Projected proportion of IPD covered by PCVs",
                            variable_name = sample_estimates$vaccine, n=NA, N=NA, sample_mean = round(sample_estimates$med*100, 1),
                            sample_confint_low = round(sample_estimates$low*100, 1), sample_confint_high = round(sample_estimates$high*100, 1),
                            pop_mean = round(population_estimates_all$med*100, 1), pop_confint_low = round(population_estimates_all$low*100, 1),
                            pop_confint_high = round(population_estimates_all$high*100, 1), popu5_mean = round(population_estimates_u5$med*100, 1),
                            popu5_confint_low = round(population_estimates_u5$low*100, 1), popu5_confint_high = round(population_estimates_u5$high*100, 1),
                            type="percentage")

table1_ipd_pcv[, n_format := ifelse(is.na(n), N, sprintf("%s/%s", n, N))]
table1_ipd_pcv[!is.na(sample_confint_low), sample_confint :=
                                        sprintf("%s - %s (IQR)", pmax(0, sample_confint_low), sample_confint_high)]
table1_ipd_pcv[, pop_confint :=
                                        sprintf("%s - %s", pmax(0, pop_confint_low), pop_confint_high)]
table1_ipd_pcv[, popu5_confint :=
                                        sprintf("%s - %s", pmax(0, popu5_confint_low), popu5_confint_high)]
table1_ipd_pcv[, sample_mean_format := ifelse(type == "percentage", sprintf("%s%%", sample_mean), sample_mean)]
table1_ipd_pcv[, pop_mean_format := ifelse(type == "percentage", sprintf("%s%%", pop_mean), pop_mean)]
table1_ipd_pcv[, popu5_mean_format := ifelse(type == "percentage", sprintf("%s%%", popu5_mean), popu5_mean)]

if(make_table) table1_ipd_pcv[, c("heading", "variable_name", "n_format",
                                  "sample_mean_format", "sample_confint",
                                  "pop_mean_format", "pop_confint",
                                  "popu5_mean_format", "popu5_confint")] %>%
  kblOut(booktabs=T, align=c("l","l","l","r","l","r","l","r","l"), linesep = "",
         col.names=c("", "", "total", "sample value", "", "population estimate", "", "population estimate (<5y)", ""),
         other_functions = list(function(x) collapse_rows(x, 1:2, row_group_label_position = "stack",
                                                          row_group_label_fonts = list(list(bold=T, italic=T),
                                                                                       list(bold=T, italic=T))),
                                function(x) kable_styling(x, latex_options = "scale_down")),
         out_name = "table1b_ipd_covered_pcv")

#############################################################################################################################
#' Section 8. Generate Supplemental Table D1: Proportion of current IPD covered by PCVs
#############################################################################################################################
table_ipd_pcv = population_estimates_byage %>%
  rbind(population_estimates_all %>% .[, agegrp := "All"]) %>%
  .[, agegrp := factor(agegrp, c("All", age_groups_sample$name))] %>%
  setorder(agegrp, vaccine) %>%
  .[, report := sprintf("%s%% (%s - %s)", round(med*100, 1), round(low*100, 1), round(high*100, 1))] %>%
  dcast(agegrp~vaccine, value.var="report")

if(make_table) table_ipd_pcv %>%
  kblOut(booktabs=T, align=c("l","l","l","l","l","l"), linesep = "",
         col.names=c("Age group", "PNEUMOSIL", "Synflorix", "Prevenar 13", "Vaxneuvance", "Prevenar 20"),
         other_functions = list(function(x) kable_styling(x, latex_options = "scale_down")),
         out_name = "tableSD1_ipd_covered_pcv")
