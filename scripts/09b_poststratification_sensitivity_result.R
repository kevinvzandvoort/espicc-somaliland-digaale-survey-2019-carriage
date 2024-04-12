#' In which output directory to store this output?
OUTPUT_DIR = "ps_age_sex"

#############################################################################################################################
#' Section 1. Generate Supplemental Table C5. Pneumococcal prevalence estimates by different post-stratification weights.
#############################################################################################################################
table_data = postratification_sensitivity_data %>%
  lapply("[[", "table_data") %>% rbindlist

#' Combine al sections in the table
table_data = table_data %>%
  .[, c("heading", "variable_name", "n", "N", "sample_mean", "pop_mean", "pop_confint_low", "pop_confint_high",
        "popu5_mean", "popu5_confint_low", "popu5_confint_high", "type")]
table_data[, n_format := ifelse(is.na(n), N, sprintf("%s/%s", n, N))]
table_data[, pop_confint := sprintf("%s - %s", pmax(0, pop_confint_low), pop_confint_high)]
table_data[, popu5_confint := sprintf("%s - %s", pmax(0, popu5_confint_low), popu5_confint_high)]
table_data[, sample_mean_format := ifelse(type == "percentage", sprintf("%s%%", sample_mean), sample_mean)]
table_data[, pop_mean_format := ifelse(type == "percentage", sprintf("%s%%", pop_mean), pop_mean)]
table_data[, popu5_mean_format := ifelse(type == "percentage", sprintf("%s%%", popu5_mean), popu5_mean)]

#' Write the table to the output folder
if(make_table) table_data[, c("heading", "variable_name", "n_format", "sample_mean_format", "pop_mean_format", "pop_confint",
                              "popu5_mean_format", "popu5_confint")] %>%
  kblOut(booktabs=T,
         align=c("l","l","l","r","r","l","r","l"),
         linesep = "",
         col.names=c("", "", "total", "sample value", "population estimate", "", "population estimate (<5y)", ""),
         other_functions = list(function(x) kable_styling(x, latex_options = "scale_down")),
         out_name = "tableSC5_post_strat_sens")

#############################################################################################################################
#' Section 2. Generate Supplemental Figure C4. Prevalence and serotype distribution by age using different weights.
#############################################################################################################################
plot_data_prevalence_overall = postratification_sensitivity_data %>%
  lapply(function(x) x[["plot_data"]][["prevalence_data"]]) %>% rbindlist

plot_data_prevalence_bytype = postratification_sensitivity_data %>%
  lapply(function(x) x[["plot_data"]][["prevalence_data_bytype"]]) %>% rbindlist

figureSC4_prevalence_byage_bystratum = plot_data_prevalence_bytype %>%
ggplot(aes(x=agegrp, y=mean, fill=factor(variable, c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT", "NVT", "NVT + NT", "NT"),
                                         c("VT", "VT + NESp", "VT + NVT", "VT + NVT + NESp", "NVT", "NVT + NESp", "NESp"))))+
  facet_wrap(nrow = 2, ncol=2, factor(stratum, c("age_sex", "age_householdsize", "age", "none"),
                    c("I", "II", "III", "IV"))~.)+
  geom_col(width=0.5)+
  scale_fill_manual(values = c("VT" = lshtm_colours$lightgreen,
                               "VT + NESp" = lshtm_colours$darkgreen,
                               "VT + NVT" = lshtm_colours$blue,
                               "VT + NVT + NESp" = lshtm_colours$darkblue,
                               "NVT" = lshtm_colours$pink,
                               "NVT + NESp" = lshtm_colours$purple,
                               "NESp" = lshtm_colours$darkgrey))+
  scale_y_continuous(labels = scales::percent)+
  theme_bw()+
  coord_cartesian(ylim=c(0,1))+
  labs(x="Age group", y="Prevalence", fill="Serotype(s) carried")+
  geom_errorbar(data = plot_data_prevalence_overall,
                aes(ymin=ci95_low, ymax=ci95_high, fill=NULL, y=NULL), width=0.05, colour=lshtm_colours$black, alpha=0.5)+
  guides(fill=guide_legend(nrow=2, byrow=FALSE))+
  theme_minimal()+
  theme_multiplot+
  theme(legend.position = "bottom")

for(ext in c("png", "pdf", "tiff", "eps"))
  ggsave(sprintf("%s/output/%s/figures/%s/figureSC4_prevalence_byage_bystratum.%s", analysis_dir, OUTPUT_DIR, ext, ext),
         plot = figureSC4_prevalence_byage_bystratum, width = plot_double_col, height = plot_double_col*0.55, units = "in", dpi = 300)