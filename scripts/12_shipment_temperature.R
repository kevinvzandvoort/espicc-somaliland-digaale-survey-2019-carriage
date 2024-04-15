#############################################################################################################################
#' Section 1. Generate Supplemental Table B1. Time until temperature exceedance.
#############################################################################################################################
data = fread("./data/other_data/shipment_temperature_data_s0_test.csv")

data_exceedance = data.table(temperature_exceeded = c(-20, -15, -10, -5, 0, 5, 10, 15),
                             hours = round(c(-20, -15, -10, -5, 0, 5, 10, 15) %>%
                                             sapply(function(C) first(data[hours_passed > 4 & y >= C, hours_passed]))))
data_exceedance[, days := round(hours/24, 1), by=.I] %>% .[]

if(make_table) data_exceedance %>%
  kblOut(booktabs=T, align=c("l","l"), linesep = "",
         col.names=c("Temperature exceeded", "Hours passed", "Days passed"),
         other_functions = list(function(x) kable_styling(x, latex_options = "scale_down")),
         out_name = "tableSB1_temperature_shipment")

#############################################################################################################################
#' Section 2. Generate Supplemental Figure B1. Temperature of test shipment over time.
#############################################################################################################################

figure_temperature_test = data %>%
  ggplot(aes(x=hours_passed, y=y))+
  geom_line(linewidth=1)+
  theme_bw()+
  scale_x_continuous(breaks = seq(0, 24*12, 24), limits = c(0, 24*12))+
  scale_y_continuous(breaks = seq(-80, 30, 5))+
  geom_hline(yintercept=0)+
  geom_hline(yintercept=c(-20), linetype=2)+
  geom_vline(xintercept=11*24, linetype=2)+
  labs(x="Hours passed", y="Temperature (°C)")+
  theme_multiplot+
  theme(axis.text.x = element_text(angle=45, hjust=1))

for(ext in c("png", "pdf", "tiff", "eps"))
  ggsave(sprintf("%s/output/%s/figures/%s/figureSB1_temperature_shipment_test.%s", analysis_dir, OUTPUT_DIR, ext, ext),
         plot = figure_temperature_test, width = plot_double_col, height = plot_single_col, units = "in", dpi = 300)

#############################################################################################################################
#' Section 3. Supplemental Figure B2. Temperature of pilot shipment over time.
#############################################################################################################################

data = fread("./data/other_data/shipment_temperature_data_s1_pilot.csv")

figure_temperature_pilot = data %>%
  ggplot(aes(x=hours_passed, y=y))+
  geom_line(size=1)+
  theme_bw()+
  scale_x_continuous(breaks = seq(0, 24*31, 24), limits = c(0, 24*31))+
  coord_cartesian(xlim=c(0, 24*29))+
  scale_y_continuous(breaks = seq(-90, 30, 5))+
  geom_hline(yintercept=0)+
  geom_hline(yintercept=c(-20), linetype=2)+
  labs(x="Hours passed", y="Temperature (°C)")+
  theme_multiplot+
  theme(axis.text.x = element_text(angle=45, hjust=1))

for(ext in c("png", "pdf", "tiff", "eps"))
  ggsave(sprintf("%s/output/%s/figures/%s/figureSB2_temperature_shipment_pilot.%s", analysis_dir, OUTPUT_DIR, ext, ext),
         plot = figure_temperature_pilot, width = plot_double_col, height = plot_single_col, units = "in", dpi = 300)
