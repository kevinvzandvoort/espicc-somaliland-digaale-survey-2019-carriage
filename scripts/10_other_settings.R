setwd(analysis_dir)

#############################################################################################################################
#' Section 1. Read data from other settings
#' - The Gambia: https://doi.org/10.1086/506941
#' - Kenya: https://doi.org/10.1016/S2214-109X(14)70224-4
#' - Malawi: https://doi.org/10.1016/j.vaccine.2018.10.021
#' - Uganda: https://doi.org/10.1016%2Fj.vaccine.2017.07.081
#############################################################################################################################
data_gambia = fread("./data/other_data/carriage_gambia.csv")
data_gambia %<>%
  .[, c("VT", "NVT", "S") := .(VT*N, NVT*N, S*N)] %>%
  melt(measure.vars=c("S", "VT", "NVT"), variable.name = "compartment", value.name = "observed") %>%
  .[, age_group := age] %>% .[, -"age"] %>%
  merge(data.table(age_group = c("<1", "1-4", "5-14", "15-39", "40+"),
                   age_low = c(0, 1, 5, 15, 40), age_high = c(1, 5, 15, 40, 60))) %>%
  dcast(...~compartment, value.var = "observed") %>% .[, any := VT+NVT] %>%
  melt(measure.vars = c("S", "VT", "NVT", "any"), variable.name="compartment", value.name="observed") %>%
  .[, c("mid", "low", "high") := binom.confint(observed, N, methods="exact")[, c("mean", "lower", "upper")]] %>%
  .[, c("compartment", "age_group", "age_low", "age_high", "mid", "low", "high")] %>% .[, country := "Gambia"]

data_kenya = fread("./data/other_data/carriage_kenya_kilifi.csv")
data_kenya %<>%
  .[, S := N - VT - NVT] %>%
  melt(measure.vars=c("S", "VT", "NVT"), variable.name = "compartment", value.name = "observed") %>%
  merge(data.table(age_group = c("<1", "1-5", "6-14", "15-20", "21-49", "50+"),
                   age_low = c(0, 1, 6, 15, 21, 50), age_high = c(1, 6, 15, 21, 50, 60)), by = "age_group") %>%
  dcast(...~compartment, value.var = "observed") %>% .[, any := VT+NVT] %>%
  melt(measure.vars = c("S", "VT", "NVT", "any"), variable.name="compartment", value.name="observed") %>%
  .[, c("mid", "low", "high") := binom.confint(observed, N, methods="exact")[, c("mean", "lower", "upper")]] %>%
  .[, c("compartment", "age_group", "age_low", "age_high", "mid", "low", "high")] %>% .[, country := "Kenya"]

data_malawi = fread("./data/other_data/carriage_malawi.csv")
data_malawi %<>%
  .[, S := N - VT - NVT] %>%
  melt(measure.vars=c("S", "VT", "NVT"), variable.name = "compartment", value.name = "observed") %>%
  .[, age_group := age]  %>%
  merge(data.table(age_group = c("<1", "1-4", "5-15", "20-35"),
                   age_low = c(0, 1, 5, 20), age_high = c(1, 5, 15, 35)), by = "age_group") %>%
  dcast(...~compartment, value.var = "observed") %>% .[, any := VT+NVT] %>%
  melt(measure.vars = c("S", "VT", "NVT", "any"), variable.name="compartment", value.name="observed") %>%
  .[, c("mid", "low", "high") := binom.confint(observed, N, methods="exact")[, c("mean", "lower", "upper")]] %>%
  .[, c("compartment", "age_group", "age_low", "age_high", "mid", "low", "high")] %>% .[, country := "Malawi"]

data_uganda = fread("./data/other_data/carriage_uganda.csv")
data_uganda %<>%
  .[, c("VT", "NVT", "S") := .(VT13, TOTALST-VT10, N-TOTALST)] %>%
  melt(measure.vars=c("S", "VT", "NVT"), variable.name = "compartment", value.name = "observed") %>%
  .[, age_group := age] %>% .[, -"age"] %>%
  merge(data.table(age_group = c("<2", "2-4", "5-14", "15+"),
                   age_low = c(0, 2, 5, 15), age_high = c(1, 5, 15, 50)), by = "age_group") %>%
  dcast(...~compartment, value.var = "observed") %>% .[, any := VT+NVT] %>%
  melt(measure.vars = c("S", "VT", "NVT", "any"), variable.name="compartment", value.name="observed") %>%
  .[, c("mid", "low", "high") := binom.confint(observed, N, methods="exact")[, c("mean", "lower", "upper")]] %>%
  .[, c("compartment", "age_group", "age_low", "age_high", "mid", "low", "high")] %>% .[, country := "Uganda"]

data_digaale = svyMean2(~dom_pcv10_type, participant_data_design_ps, na.rm=T, by="agegrp") %>%
  .[option != "NT"] %>%
  .[, c("agegrp", "option", "mean", "ci95_low", "ci95_high")] %>%
  .[, c("compartment", "age_group", "mid", "low", "high", "country") :=
      .(option, agegrp, mean, ci95_low, ci95_high, "Digaale IDP camp (SLD)")] %>% 
  .[, low := pmax(0, low)] %>%
  .[, c("compartment", "age_group", "mid", "low", "high", "country")] %>%
  rbind(svyMean2(~pneu_carr_final, participant_data_design_ps, na.rm=T, by="agegrp") %>%
          .[, option := "any"] %>%
          .[, c("agegrp", "option", "mean", "ci95_low", "ci95_high")] %>%
          .[, c("compartment", "age_group", "mid", "low", "high", "country") :=
              .(option, agegrp, mean, ci95_low, ci95_high, "Digaale IDP camp (SLD)")] %>% 
          .[, low := pmax(0, low)] %>%
          .[, c("compartment", "age_group", "mid", "low", "high", "country")]) %>%
  merge(age_groups_sample %>% copy %>% setNames(c("age_low", "age_high", "age_group")), by="age_group", all.x=TRUE) %>%
  .[, age_high := pmin(60, age_high)]
data_digaale

#############################################################################################################################
#' Section 2. Generate Supplemental Figure C5. Prevalence by age compared to different settings.
#############################################################################################################################  
data_carriage = rbind(data_gambia, data_kenya, data_malawi, data_uganda, data_digaale)

(figureSC5_prevalence_differentsettings = data_carriage[compartment != "S"] %>%
  melt(measure.vars = c("mid", "low", "high")) %>%
  dcast(... ~ variable, value.var="value") %>%
  .[, setting := ifelse(country == "Digaale IDP camp (SLD)", "Digaale IDP camp (SLD)", "Other settings")] %>%
  rbind(copy(.) %>% .[, setting := "Combined"]) %>%
  .[setting == "Combined"] %>%
  ggplot(aes(x = age_low, xmin=age_low, xend = age_high, xmax=age_high,
             ymin = low, ymax = high, fill = country, colour = country))+
  facet_grid(. ~ factor(compartment, c("any", "VT", "NVT"),
                    c("Any pneumococci carrier", "Vaccine Type carrier", "Non-Vaccine Type carrier")))+
  geom_rect(alpha=0.2, colour=NA)+
  geom_segment(aes(y=mid, yend=mid), linewidth=1)+
  labs(x="Age (years)", y="Prevalence", fill="Country carriage data", colour="Country carriage data")+
  scale_y_continuous(labels=scales::percent, limits = c(0, 1))+
  scale_x_continuous(breaks = seq(0, 60, 10))+
  scale_colour_manual(values = c("Kenya" = lshtm_colours$lightgreen, "Gambia" = lshtm_colours$pink,
                                 "Malawi" = lshtm_colours$blue, "Uganda" = lshtm_colours$yellow,
                                 "Digaale IDP camp (SLD)" = lshtm_colours$black))+
  scale_fill_manual(values = c("Kenya" = lshtm_colours$lightgreen, "Gambia" = lshtm_colours$pink,
                               "Malawi" = lshtm_colours$blue, "Uganda" = lshtm_colours$yellow,
                               "Digaale IDP camp (SLD)" = lshtm_colours$darkgrey))+
  theme_minimal()+
  theme_multiplot)

for(ext in c("png", "pdf", "tiff", "eps"))
  ggsave(sprintf("%s/output/%s/figures/%s/figureSC5_prevalence_other_settings.%s", analysis_dir, OUTPUT_DIR, ext, ext),
         plot = figureSC5_prevalence_differentsettings, width = plot_double_col*1.2, height = plot_single_col, units = "in", dpi = 300)