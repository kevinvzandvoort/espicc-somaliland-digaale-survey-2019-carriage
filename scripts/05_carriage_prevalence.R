#############################################################################################################################
#' Section 1. Carriage prevalence by age - VTs defined as Pneumosil VT
#############################################################################################################################
plot_byage_pneumosil = svyMean2(~pneumosil_type, participant_data_design_ps, na.rm=T, by="agegrp") %>%
  .[, c("agegrp", "option", "mean")] %>%
  .[option != "S"] %>%
  .[, variable := option] %>%
  .[] %>%
  ggplot(aes(x=agegrp, y=mean, fill=factor(variable, c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT", "NVT", "NVT + NT", "NT"),
                                           c("VT", "VT + NESp", "VT + NVT", "VT + NVT + NESp", "NVT", "NVT + NESp", "NESp"))))+
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
  geom_errorbar(data = svyMean2(~pneu_carr_final, participant_data_design_ps, na.rm=T, by="agegrp"),
                aes(ymin=ci95_low, ymax=ci95_high, fill=NULL, y=NULL), width=0.05, colour=lshtm_colours$black, alpha=0.5)+
  guides(fill=guide_legend(nrow=2, byrow=FALSE))+
  theme_minimal()+
  theme_multiplot+
  theme(legend.position = "bottom")

#############################################################################################################################
#' Section 2.  Generate Figure 3: Prevalence and serotype distribution by age.
#############################################################################################################################
figure3_prevalence_byage_pneumosil = plot_byage_pneumosil
for(ext in c("png", "pdf", "tiff", "eps"))
  ggsave(sprintf("%s/output/%s/figures/%s/figure3_prevalence_byage_pneumosil.%s", analysis_dir, OUTPUT_DIR, ext, ext),
         plot = figure3_prevalence_byage_pneumosil, width = plot_single_col*1.75, height = plot_single_col, units = "in", dpi = 300)

#############################################################################################################################
#' Section 3.  Prevalence by age age, using all PCVs.
#############################################################################################################################
figureSC3_prevalence_byage =
  svyMean2(~pneumosil_type, participant_data_design_ps, na.rm=T, by="agegrp") %>%
  .[, c("agegrp", "option", "mean")] %>%
  .[option != "S"] %>%
  .[, c("variable", "type", "type2") := .(option, "PNEUMOSIL", "Multiple carriage")] %>%
  rbind(svyMean2(~pcv10_type, participant_data_design_ps, na.rm=T, by="agegrp") %>%
          .[, c("agegrp", "option", "mean")] %>%
          .[option != "S"] %>%
          .[, c("variable", "type", "type2") := .(option, "PCV10", "Multiple carriage")]) %>%
  rbind(svyMean2(~pcv13_type, participant_data_design_ps, na.rm=T, by="agegrp") %>%
          .[, c("agegrp", "option", "mean")] %>%
          .[option != "S"] %>%
          .[, c("variable", "type", "type2") := .(option, "PCV13", "Multiple carriage")]) %>%
  rbind(svyMean2(~pcv15_type, participant_data_design_ps, na.rm=T, by="agegrp") %>%
          .[, c("agegrp", "option", "mean")] %>%
          .[option != "S"] %>%
          .[, c("variable", "type", "type2") := .(option, "PCV15", "Multiple carriage")]) %>%
  rbind(svyMean2(~pcv20_type, participant_data_design_ps, na.rm=T, by="agegrp") %>%
          .[, c("agegrp", "option", "mean")] %>%
          .[option != "S"] %>%
          .[, c("variable", "type", "type2") := .(option, "PCV20", "Multiple carriage")]) %>%
  rbind(svyMean2(~dom_pneumosil_type, participant_data_design_ps, na.rm=T, by="agegrp") %>%
          .[, c("agegrp", "option", "mean")] %>%
          .[option != "S"] %>%
          .[, c("variable", "type", "type2") := .(option, "PNEUMOSIL", "Dominant carriage")]) %>%
  rbind(svyMean2(~dom_pcv10_type, participant_data_design_ps, na.rm=T, by="agegrp") %>%
          .[, c("agegrp", "option", "mean")] %>%
          .[option != "S"] %>%
          .[, c("variable", "type", "type2") := .(option, "PCV10", "Dominant carriage")]) %>%
  rbind(svyMean2(~dom_pcv13_type, participant_data_design_ps, na.rm=T, by="agegrp") %>%
          .[, c("agegrp", "option", "mean")] %>%
          .[option != "S"] %>%
          .[, c("variable", "type", "type2") := .(option, "PCV13", "Dominant carriage")]) %>%
  rbind(svyMean2(~dom_pcv15_type, participant_data_design_ps, na.rm=T, by="agegrp") %>%
          .[, c("agegrp", "option", "mean")] %>%
          .[option != "S"] %>%
          .[, c("variable", "type", "type2") := .(option, "PCV15", "Dominant carriage")]) %>%
  rbind(svyMean2(~dom_pcv20_type, participant_data_design_ps, na.rm=T, by="agegrp") %>%
          .[, c("agegrp", "option", "mean")] %>%
          .[option != "S"] %>%
          .[, c("variable", "type", "type2") := .(option, "PCV20", "Dominant carriage")]) %>%
  ggplot(aes(x=agegrp, y=mean, fill=factor(variable, c("VT", "VT + NT", "VT + NVT", "VT + NVT + NT", "NVT", "NVT + NT", "NT"),
                                           c("VT", "VT + NESp", "VT + NVT", "VT + NVT + NESp", "NVT", "NVT + NESp", "NESp"))))+
  facet_grid(factor(type, c("PNEUMOSIL", "PCV10", "PCV13", "PCV15", "PCV20"),
                    c("PNEUMOSIL", "Synflorix", "Prevenar13", "Vaxneuvance", "Prevenar20"))~type2)+
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
  geom_errorbar(data = svyMean2(~pneu_carr_final, participant_data_design_ps, na.rm=T, by="agegrp"),
                aes(ymin=ci95_low, ymax=ci95_high, fill=NULL, y=NULL), width=0.05, colour=lshtm_colours$black, alpha=0.5)+
  guides(fill=guide_legend(nrow=2, byrow=FALSE))+
  theme_minimal()+
  theme_multiplot+
  theme(legend.position = "bottom")

#############################################################################################################################
#' Section 4.  Generate Supplemental Figure C3. Prevalence and serotype distribution by age.
#############################################################################################################################
for(ext in c("png", "pdf", "tiff", "eps"))
  ggsave(sprintf("%s/output/%s/figures/%s/figureSC3_prevalence_byage.%s", analysis_dir, OUTPUT_DIR, ext, ext),
         plot = figureSC3_prevalence_byage, width = plot_double_col, height = plot_double_col, units = "in", dpi = 300)
