#############################################################################################################################
#' Section 1. Create data for serotype distribution plots, showing frequency of serotypes identified in NP swabs
#############################################################################################################################
plot_data = lab_data %>%
  melt(measure.vars = list(serotype = sprintf("sero%s_final", 1:5),
                           relative_abundance = sprintf("ra_sero%s_final", 1:5),
                           absolute_abundance = sprintf("dens_sero%s", 1:5))) %>%
  .[!is.na(serotype) & serotype != "", .(ra = sum(relative_abundance), N = .N, dens=sum(absolute_abundance),
                        domN = sum(variable == 1)), by=serotype] %>%
  .[, c("mid_ra", "low_ra", "high_ra") := binom.confint(ra, sum(ra), methods="exact") %>%
      .[, c("mean", "lower", "upper")]] %>%
  .[, c("mid_abs", "low_abs", "high_abs") := binom.confint(N, sum(N), methods="exact") %>%
      .[, c("mean", "lower", "upper")]] %>%
  .[, c("mid_dom", "low_dom", "high_dom") := binom.confint(domN, sum(domN), methods="exact") %>% 
      .[, c("mean", "lower", "upper")]] %>%
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
                                                    ifelse(pcv20, "pcv20", "NVT"))))))] %>%
  .[, vaccine := ifelse(grepl("NT", serotype, fixed=TRUE), "NT", vaccine)] %>%
  .[, c("mid", "low", "high", "type") := .(mid_ra, low_ra, high_ra, "ra")] %>% .[, -c("mid_ra", "low_ra", "high_ra")] %>%
  rbind(copy(.) %>% .[, c("mid", "low", "high", "type") := .(mid_abs, low_abs, high_abs, "absolute")] %>%
          rbind(copy(.) %>% .[, c("mid", "low", "high", "type") := .(mid_dom, low_dom, high_dom, "dominant")])) %>%
  .[, -c("mid_abs", "low_abs", "high_abs", "mid_dom", "low_dom", "high_dom")] %>%
  .[, total := ifelse(type == "ra", sum(ra), ifelse(type == "absolute", sum(N),
                                                    ifelse(type == "dominant", sum(domN), sum(dens)))), by="type"]

#############################################################################################################################
#' Section 2. Plot serotype distribution, weighted by relative abundance
#############################################################################################################################
plot_distribution_ra = plot_data[type == "ra"] %>%
  ggplot(mapping = aes(x=factor(serotype, rev(plot_data[type == "ra"] %>% setorder(mid) %>% .[, serotype])),
                       fill=factor(vaccine, c("pcv20_15_13_10_pneumosil", "pcv20_15_13_pneumosil",
                                              "pcv20_15_13_10", "pcv20_15_13", "pcv20_15", "pcv20", "NVT", "NT"),
                                   c("All", "Prevenar20, Vaxneuvance, Prevenar13, PNEUMOSIL",
                                     "Prevenar20, Vaxneuvance, Prevenar13, Synflorix",
                                     "Prevenar20, Vaxneuvance, Prevenar13", "Prevenar20, Vaxneuvance",
                                     "Prevenar20", "NVT", "NESp")),
                       y=mid, ymin=low, ymax=high))+
  geom_col()+
  geom_errorbar(width=0.15, colour=lshtm_colours$black, alpha=0.5)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "bottom")+
  labs(x="Serotype", y="Prevalence", fill="Serotype(s)")+
  scale_y_continuous(labels = scales::percent)+
  scale_fill_manual(values = c("All" = lshtm_colours$lightgreen,
                               "Prevenar20, Vaxneuvance, Prevenar13, PNEUMOSIL" = lshtm_colours$blue,
                               "Prevenar20, Vaxneuvance, Prevenar13, Synflorix" = lshtm_colours$purple,
                               "Prevenar20, Vaxneuvance, Prevenar13" = lshtm_colours$pink,
                               #"Prevenar20, Vaxneuvance" = lshtm_colours$pink,
                               "Prevenar20" = lshtm_colours$yellow,
                               "NVT" = lshtm_colours$lightgrey,
                               "NESp" = lshtm_colours$darkgrey))+
  guides(fill=guide_legend(nrow=5, ncol=2, byrow=FALSE))+
  theme_minimal()+
  theme_multiplot+
  theme(legend.position = "bottom", axis.text.x = element_text(angle=45, hjust=1))

#############################################################################################################################
#' Section 3. Plot serotype distribution for dominant serotypes only
#############################################################################################################################
plot_distribution_dominant = plot_data[type == "dominant"] %>%
  ggplot(mapping = aes(x=factor(serotype, rev(plot_data[type == "ra"] %>% setorder(mid) %>% .[, serotype])),
                       fill=factor(vaccine, c("pcv20_15_13_10_pneumosil", "pcv20_15_13_pneumosil",
                                              "pcv20_15_13_10", "pcv20_15_13", "pcv20_15", "pcv20", "NVT", "NT"),
                                   c("All", "Prevenar20, Vaxneuvance, Prevenar13, PNEUMOSIL",
                                     "Prevenar20, Vaxneuvance, Prevenar13, Synflorix",
                                     "Prevenar20, Vaxneuvance, Prevenar13", "Prevenar20, Vaxneuvance",
                                     "Prevenar20", "NVT", "NESp")),
                       y=mid, ymin=low, ymax=high))+
  geom_col()+
  geom_errorbar(width=0.15, colour=lshtm_colours$black, alpha=0.5)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "bottom")+
  labs(x="Serotype", y="Prevalence", fill="Serotype(s)")+
  scale_y_continuous(labels = scales::percent,
                     sec.axis = sec_axis(trans = function(x) x * plot_data[type == "dominant", sum(domN)]))+
  scale_fill_manual(values = c("All" = lshtm_colours$lightgreen,
                               "Prevenar20, Vaxneuvance, Prevenar13, PNEUMOSIL" = lshtm_colours$blue,
                               "Prevenar20, Vaxneuvance, Prevenar13, Synflorix" = lshtm_colours$purple,
                               "Prevenar20, Vaxneuvance, Prevenar13" = lshtm_colours$pink,
                               #"Prevenar20, Vaxneuvance" = lshtm_colours$pink,
                               "Prevenar20" = lshtm_colours$yellow,
                               "NVT" = lshtm_colours$lightgrey,
                               "NESp" = lshtm_colours$darkgrey))+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))+
  theme_minimal()+
  theme_multiplot+
  theme(legend.position = "bottom", axis.text.x = element_text(angle=45, hjust=1))

#############################################################################################################################
#' Section 4. Plot unweighted serotype distribution
#############################################################################################################################
plot_distribution_N = plot_data[type == "absolute"] %>%
  ggplot(mapping = aes(x=factor(serotype, rev(plot_data[type == "ra"] %>% setorder(mid) %>% .[, serotype])),
                       fill=factor(vaccine, c("pcv20_15_13_10_pneumosil", "pcv20_15_13_pneumosil",
                                              "pcv20_15_13_10", "pcv20_15_13", "pcv20_15", "pcv20", "NVT", "NT"),
                                   c("All", "Prevenar20, Vaxneuvance, Prevenar13, PNEUMOSIL",
                                     "Prevenar20, Vaxneuvance, Prevenar13, Synflorix",
                                     "Prevenar20, Vaxneuvance, Prevenar13", "Prevenar20, Vaxneuvance",
                                     "Prevenar20", "NVT", "NESp")),
                       y=mid, ymin=low, ymax=high))+
  geom_col()+
  geom_errorbar(width=0.15, colour=lshtm_colours$black, alpha=0.5)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "bottom")+
  labs(x="Serotype", y="Prevalence", fill="Serotype(s)")+
  scale_y_continuous(labels = scales::percent, sec.axis = sec_axis(trans = function(x) x * plot_data[type == "absolute", sum(domN)]))+
  scale_fill_manual(values = c("All" = lshtm_colours$lightgreen,
                               "Prevenar20, Vaxneuvance, Prevenar13, PNEUMOSIL" = lshtm_colours$blue,
                               "Prevenar20, Vaxneuvance, Prevenar13, Synflorix" = lshtm_colours$purple,
                               "Prevenar20, Vaxneuvance, Prevenar13" = lshtm_colours$pink,
                               #"Prevenar20, Vaxneuvance" = lshtm_colours$pink,
                               "Prevenar20" = lshtm_colours$yellow,
                               "NVT" = lshtm_colours$lightgrey,
                               "NESp" = lshtm_colours$darkgrey))+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))+
  theme_minimal()+
  theme_multiplot+
  theme(legend.position = "bottom", axis.text.x = element_text(angle=45, hjust=1))

#############################################################################################################################
#' Section 5. generates Figure 2: Pneumococcal serotype distribution.
#############################################################################################################################
figure1_serotype_distribution_ra = plot_distribution_ra
for(ext in c("png", "pdf", "tiff", "eps"))
  ggsave(sprintf("%s/output/%s/figures/%s/figure2_serotype_distribution_ra.%s", analysis_dir, OUTPUT_DIR, ext, ext),
         plot = figure1_serotype_distribution_ra, width = plot_double_col, height = plot_double_col/1.6, units = "in", dpi = 300)

#############################################################################################################################
#' Section 6. generates Supplemental Figure C1. Pneumococcal serotype distribution (alternative distributions)
#############################################################################################################################
figureSC1_serotype_distribution = (plot_distribution_ra+ggtitle(label = "All serotypes - weighted by relative abundance"))/
  (plot_distribution_dominant+ggtitle(label = "Dominant serotypes"))/
  (plot_distribution_N+ggtitle(label = "All serotypes"))+
  plot_layout(guides = "collect")+
  plot_annotation(tag_levels = 'A')&
  guides(fill=guide_legend(nrow=5, ncol=2, byrow=FALSE))&
  theme(legend.position = "bottom")
for(ext in c("png", "pdf", "tiff", "eps"))
  ggsave(sprintf("%s/output/%s/figures/%s/figureSC1_serotype_distribution.%s", analysis_dir, OUTPUT_DIR, ext, ext),
         plot = figureSC1_serotype_distribution, width = plot_double_col, height = (plot_double_col), units = "in", dpi = 300)

#############################################################################################################################
#' Section 6. generates Supplemental Figure C2. Pneumococcal serotype distribution (by age).
#############################################################################################################################
plot_data_byage = lab_data %>%
  merge(participant_data[, c("pid", "participant_age_y")], by="pid") %>%
  .[!is.na(participant_age_y)] %>%
  .[, u1 := participant_age_y < 1] %>%
  .[, u5 := participant_age_y < 5] %>%
  melt(measure.vars = list(serotype = sprintf("sero%s_final", 1:5),
                           relative_abundance = sprintf("ra_sero%s_final", 1:5),
                           absolute_abundance = sprintf("dens_sero%s", 1:5))) %>%
  .[!is.na(serotype) & serotype != "", .(ra = sum(relative_abundance), N = .N, dens=sum(absolute_abundance),
                        domN = sum(variable == 1)), by=c("serotype", "u5")] %>%
  .[, c("mid_ra", "low_ra", "high_ra") := binom.confint(ra, sum(ra), methods="exact") %>%
      .[, c("mean", "lower", "upper")]] %>%
  .[, c("mid_abs", "low_abs", "high_abs") := binom.confint(N, sum(N), methods="exact") %>%
      .[, c("mean", "lower", "upper")]] %>%
  .[, c("mid_dom", "low_dom", "high_dom") := binom.confint(domN, sum(domN), methods="exact") %>% 
      .[, c("mean", "lower", "upper")]] %>%
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
                                                    ifelse(pcv20, "pcv20", "NVT"))))))] %>%
  .[, vaccine := ifelse(grepl("NT", serotype, fixed=TRUE), "NT", vaccine)] %>%
  .[, c("mid", "low", "high", "type") := .(mid_ra, low_ra, high_ra, "ra")] %>% .[, -c("mid_ra", "low_ra", "high_ra")] %>%
  rbind(copy(.) %>% .[, c("mid", "low", "high", "type") := .(mid_abs, low_abs, high_abs, "absolute")] %>%
          rbind(copy(.) %>% .[, c("mid", "low", "high", "type") := .(mid_dom, low_dom, high_dom, "dominant")])) %>%
  .[, -c("mid_abs", "low_abs", "high_abs", "mid_dom", "low_dom", "high_dom")] %>%
  .[, total := ifelse(type == "ra", sum(ra), ifelse(type == "absolute", sum(N),
                                                    ifelse(type == "dominant", sum(domN), sum(dens)))), by="type"]

figureSC2_serotype_distribution_ra_bygage = plot_data_byage %>% .[, age := ifelse(u5, "<5", "≥5")] %>%
  .[type == "ra"] %>%
  ggplot(mapping = aes(x=factor(serotype, c(rev(plot_data_byage[age == "<5" & type == "ra"] %>% setorder(mid) %>% .[, serotype]),
                                            (plot_data_byage[age != "<5", unique(serotype)] %>% subset(!. %in% plot_data_byage[age == "<5", serotype])))),
                       fill=factor(vaccine, c("pcv20_15_13_10_pneumosil", "pcv20_15_13_pneumosil",
                                              "pcv20_15_13_10", "pcv20_15_13", "pcv20_15", "pcv20", "NVT", "NT"),
                                   c("All", "Prevenar20, Vaxneuvance, Prevenar13, PNEUMOSIL",
                                     "Prevenar20, Vaxneuvance, Prevenar13, Synflorix",
                                     "Prevenar20, Vaxneuvance, Prevenar13", "Prevenar20, Vaxneuvance",
                                     "Prevenar20", "NVT", "NESp")),
                       y=mid, ymin=low, ymax=high))+
  geom_col()+
  geom_errorbar(width=0.15, colour=lshtm_colours$black, alpha=0.5)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "bottom")+
  labs(x="Serotype", y="Prevalence", fill="Serotype(s)")+
  scale_y_continuous(labels = scales::percent)+
  scale_fill_manual(values = c("All" = lshtm_colours$lightgreen,
                               "Prevenar20, Vaxneuvance, Prevenar13, PNEUMOSIL" = lshtm_colours$blue,
                               "Prevenar20, Vaxneuvance, Prevenar13, Synflorix" = lshtm_colours$purple,
                               "Prevenar20, Vaxneuvance, Prevenar13" = lshtm_colours$pink,
                               #"Prevenar20, Vaxneuvance" = lshtm_colours$pink,
                               "Prevenar20" = lshtm_colours$yellow,
                               "NVT" = lshtm_colours$lightgrey,
                               "NESp" = lshtm_colours$darkgrey))+
  guides(fill=guide_legend(nrow=5, ncol=2, byrow=FALSE))+
  facet_grid(age~., scales = "free")+
  theme_minimal()+
  theme_multiplot+
  theme(legend.position = "bottom", axis.text.x = element_text(angle=45, hjust=1))
for(ext in c("png", "pdf", "tiff", "eps"))
  ggsave(sprintf("%s/output/%s/figures/%s/figureSC2_serotype_distribution_ra_bygage.%s", analysis_dir, OUTPUT_DIR, ext, ext),
         plot = figureSC2_serotype_distribution_ra_bygage, width = plot_double_col, height = (plot_double_col*0.7), units = "in", dpi = 300)