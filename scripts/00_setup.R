#' Set global formatting options
options(knitr.kable.NA = "")
options(scipen=999)
set.seed(10)

#' Set global plot styling options
page_outer_margin = 0.523 #inch
page_inner_margin = page_outer_margin/2 #inch
page_width = 8 #inch
page_height = 11 #inch
plot_single_col = (page_width - page_outer_margin*2 - page_inner_margin)/2 #inch
plot_double_col = page_width - page_outer_margin*2 #inch

#' Additional ggplot theme options to be used when combining plots (utilizing patchwork)
theme_multiplot = ggplot2::theme(legend.title = element_text(size=9), legend.text = element_text(size = 8),
                                 legend.box.background = element_rect(fill="#FFFFFF", colour="#000000", linewidth = 0.2),
                                 legend.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "cm"),
                                 legend.key.size = unit(0.3, 'cm'), axis.title = element_text(size=10),
                                 axis.text = element_text(size = 8),
                                 plot.title = element_text(size=10, face = "bold", hjust=0.62),
                                 plot.title.position = "plot", plot.tag = element_text(size=10),
                                 plot.tag.position = c(0, 0.99), strip.text = element_text(face = "plain", size=9))

lshtm_colours = list(blue = "#00AEC7", red = "#FE5000", yellow = "#FFB81C", pink = "#FFABBA", darkgreen = "#0D5257",
                     lightgreen = "#00BF6F", darkblue = "#1E22AA", purple = "#621244", white = "#FFFFFF",
                     lightgrey = "#D9E1E2", darkgrey = "#A2ACAB", black = "#000000")

#' Create output folders
for(d in c("/output", "/output/%s", "output/%s/tables", "output/%s/tables/pdf", "output/%s/tables/html", "output/%s/figures",
           "output/%s/figures/png", "output/%s/figures/pdf", "output/%s/figures/tiff", "output/%s/figures/eps") %>%
    sprintf(paste0("ps_", tolower(POSTSTRATIFICATION_STRATA))))
  if(!dir.exists(sprintf("%s/%s", analysis_dir, d))) sprintf("%s/%s", analysis_dir, d) %>% dir.create()

OUTPUT_DIR = paste0("ps_", tolower(POSTSTRATIFICATION_STRATA))

#' Load helper functions
source(sprintf("%s/scripts/functions.R", analysis_dir))

#' Load the cleaned and anonymized data
#' - Nb some data are aggregated, to preserve anonymity of the participants
source(sprintf("%s/scripts/00a_load_data.R", analysis_dir))
#' - Pre-process data, create new variables
source(sprintf("%s/scripts/00b_process_data.R", analysis_dir))
#' - Get estimated total population size in Digaale, used for finite population corrections (FPC)
source(sprintf("%s/scripts/00c_fpc_correction.R", analysis_dir))
#' - Create Survey package objects, poststratifying the data, and applying FPCs
source(sprintf("%s/scripts/00d_poststratification.R", analysis_dir))
