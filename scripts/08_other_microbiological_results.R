#############################################################################################################################
#' Section 1. Generates Supplemental Table C1. Other microbiological results.
#############################################################################################################################
#' Multiple carriage
tableSC1.1 = rbind(data.table(variable_name = "Median number of serotypes (in all samples with pneumococci)",
                            n = NA_integer_,
                            N = lab_data %>% .[!is.na(num_sero_final) & num_sero_final >= 1, .N],
                            sample_mean = lab_data %>% .[!is.na(num_sero_final) & num_sero_final >= 1] %>% .[, num_sero_final] %>% quantile(c(0.5)),
                            sample_confint_low = lab_data %>% .[!is.na(num_sero_final) & num_sero_final >= 1] %>% .[, num_sero_final] %>% quantile(c(0.25)),
                            sample_confint_high = lab_data %>% .[!is.na(num_sero_final) & num_sero_final >= 1] %>% .[, num_sero_final] %>% quantile(c(0.75)),
                            type = "IQR"),
                 lab_data %>%
                   .[!is.na(num_sero_final)] %>%
                   .[, .(num_sero_final = sum(num_sero_final >= 2), total = sum(num_sero_final >= 1))] %>%
                   .[, binom::binom.confint(num_sero_final, total, methods = "exact")] %>% as.data.table %>%
                   .[, -"method"] %>% setNames(c("n", "N", "sample_mean", "sample_confint_low", "sample_confint_high")) %>%
                   .[, c("variable_name", "type") := .("Multiple serotype carriage", "percentage")]) %>%
  .[, heading := "Multiple serotype carriage"]

#' Median density
tableSC1.2 = rbind(data.table(variable_name = "Median pneumococcal density (total for each sample)",
                            n = NA_integer_,
                            N = lab_data %>% .[!is.na(num_sero_final) & num_sero_final > 0, .N],
                            sample_mean = lab_data %>% .[!is.na(num_sero_final) & num_sero_final > 0] %>% .[, pneu_dens_final] %>% log10 %>% quantile(c(0.5)),
                            sample_confint_low = lab_data %>% .[!is.na(num_sero_final) & num_sero_final > 0] %>% .[, pneu_dens_final] %>% log10 %>% quantile(c(0.25)),
                            sample_confint_high = lab_data %>% .[!is.na(num_sero_final) & num_sero_final > 0] %>% .[, pneu_dens_final] %>% log10 %>% quantile(c(0.75)),
                            type = "IQR"),
                 data.table(variable_name = "Median pneumococcal density (samples with only 1 serotype)",
                            n = NA_integer_,
                            N = lab_data %>% .[!is.na(num_sero_final) & num_sero_final == 1, .N],
                            sample_mean = lab_data %>% .[!is.na(num_sero_final) & num_sero_final == 1] %>% .[, pneu_dens_final] %>% log10 %>% quantile(c(0.5)),
                            sample_confint_low = lab_data %>% .[!is.na(num_sero_final) & num_sero_final == 1] %>% .[, pneu_dens_final] %>% log10 %>% quantile(c(0.25)),
                            sample_confint_high = lab_data %>% .[!is.na(num_sero_final) & num_sero_final == 1] %>% .[, pneu_dens_final] %>% log10 %>% quantile(c(0.75)),
                            type = "IQR"),
                 data.table(variable_name = "Median pneumococcal density (samples with > 1 serotype)",
                            n = NA_integer_,
                            N = lab_data %>% .[!is.na(num_sero_final) & num_sero_final > 1, .N],
                            sample_mean = lab_data %>% .[!is.na(num_sero_final) & num_sero_final > 1] %>% .[, pneu_dens_final] %>% log10 %>% quantile(c(0.5)),
                            sample_confint_low = lab_data %>% .[!is.na(num_sero_final) & num_sero_final > 1] %>% .[, pneu_dens_final] %>% log10 %>% quantile(c(0.25)),
                            sample_confint_high = lab_data %>% .[!is.na(num_sero_final) & num_sero_final > 1] %>% .[, pneu_dens_final] %>% log10 %>% quantile(c(0.75)),
                            type = "IQR"),
                 {
                   x = lab_data %>% melt(measure.vars = sprintf("dens_sero%s", 1:5)) %>% .[!is.na(value)]
                   data.table(variable_name = "Median pneumococcal density (per serotype)",
                              n = NA_integer_,
                              N = x[, .N],
                              sample_mean = x %>% .[, value] %>% log10 %>% quantile(c(0.5)),
                              sample_confint_low = x %>% .[, value] %>% log10 %>% quantile(c(0.25)),
                              sample_confint_high = x %>% .[, value] %>% log10 %>% quantile(c(0.75)),
                              type = "IQR")
                 }) %>% .[, heading := "Pneumococcal density"]

#' Non pneumococcal species
tableSC1.3 = lab_data %>%
  .[!is.na(other_spp_final)] %>%
  .[, .(other_spp = sum(other_spp_final == 1), total = .N)] %>%
  .[, binom::binom.confint(other_spp, total, methods = "exact")] %>% as.data.table %>%
  .[, -"method"] %>% setNames(c("n", "N", "sample_mean", "sample_confint_low", "sample_confint_high")) %>%
  .[, c("variable_name", "type") := .("Other species detected by microarray", "percentage")] %>%
  .[, heading := "Other species detected by microarray"]

#' Prevalence of resistant genes
tableSC1.4 = rbind(lab_data[pneu_carr_final == 1 & other_spp_final == 0 & num_sero_final == 1] %>%
                   .[, .(resistant = sum(num_abrgen_final >= 1), total = .N)] %>%
                   .[, binom::binom.confint(resistant, total, methods = "exact")] %>% as.data.table %>%
                   .[, -"method"] %>% setNames(c("n", "N", "sample_mean", "sample_confint_low", "sample_confint_high")) %>%
                   .[, c("variable_name", "type") := .("Prevalence of resistant genes", "percentage")] %>%
                   .[, subheading := ""],
                 #lab_data[pneu_carr_final == 1 & other_spp_final == 0 & num_sero_final == 1 & num_abrgen_final >= 1] %>%
                 lab_data[pneu_carr_final == 1 & other_spp_final == 0 & num_sero_final == 1] %>%
                   .[, total_resistant := .N] %>%
                   melt(measure.vars = c("apha3_final", "cat_final", "ermb_final", "ermc_final", "mefa_final", "tetk_final",
                                         "tetl_final", "tetm_final", "teto_final", "sat4_final")) %>%
                   .[, .(positive = sum(value), total = sum(!is.na(value))), by="variable"] %>%
                   .[, binom::binom.confint(positive, total, methods = "exact"), by="variable"] %>%
                   .[, gene := factor(strsplit(as.character(variable), split = "_", fixed=TRUE) %>% sapply("[[", 1),
                                      c("apha3", "cat", "ermb", "ermc", "mefa", "tetk", "tetl", "tetm", "teto", "sat4"),
                                      c("aphA3", "cat", "ermB", "ermC", "mefA", "tetK", "tetL", "tetM", "tetO", "sat4"))] %>%
                   .[, -c("variable", "method")] %>%
                   setNames(c("n", "N", "sample_mean", "sample_confint_low", "sample_confint_high", "variable_name")) %>%
                   .[, c("type") := .("percentage")] %>% .[, subheading := "Genes identified"]) %>%
  .[, heading := "Resistant genes"]

#' Combine al sections in the table
tableSC1_other_microbiological_results =
  rbind(tableSC1.1, tableSC1.2, tableSC1.3, tableSC1.4, fill=T, use.names=T) %>%
  .[, obs := ifelse(is.na(n), as.character(N), sprintf("%s/%s", n, N))] %>%
  .[, c("heading", "subheading", "variable_name", "obs", "sample_mean", "sample_confint_low", "sample_confint_high", "type")]
tableSC1_other_microbiological_results[type == "percentage", c("sample_mean", "sample_confint_low", "sample_confint_high") :=
                                         .(sample_mean*100, sample_confint_low*100, sample_confint_high*100)]
tableSC1_other_microbiological_results[, confint := sprintf("%s - %s%s", pmax(0, round(sample_confint_low, ifelse(heading == "Pneumococcal density", 2, 0))),
                                                            round(sample_confint_high, ifelse(heading == "Pneumococcal density", 2, 0)),
                                                            ifelse(type == "IQR", " (IQR)", ""))]
tableSC1_other_microbiological_results[, mean_format := ifelse(type == "percentage", sprintf("%s%%", round(sample_mean, ifelse(heading == "Pneumococcal density", 2, 0))),
                                                              round(sample_mean, ifelse(heading == "Pneumococcal density", 2, 0)))]

#' Write the table to the output folder
if(make_table) tableSC1_other_microbiological_results[, c("heading", "subheading", "variable_name", "obs", "mean_format",
                                                         "confint")] %>%
  kblOut(booktabs=T, align=c("l","l","l","l","r","l"), linesep = "", col.names=c("", "", "", "total", "value", ""),
         other_functions = list(function(x) collapse_rows(x, 1:2, row_group_label_position = "stack",
                                                          row_group_label_fonts = list(list(bold=T, italic=T),
                                                                                       list(bold=F, italic=T))),
                                function(x) kable_styling(x, latex_options = "scale_down")),
         out_name = "tableSC1_other_microbiological_results")

#############################################################################################################################
#' Section 2. Generates Supplemental Table C2. Association between serotype and dominant carriage.
#############################################################################################################################
tableSC2 = rbind(
  {
    x = lab_data %>% .[pneu_carr_final == 1] %>%
      melt(measure.vars = c("sero1_final", "sero2_final", "sero3_final", "sero4_final", "sero5_final"),
           variable.name="st_index", value.name="st") %>%
      .[st != ""] %>%
      .[, .(sttype = ifelse(st %in% serotypes_by_vaccine[pneumosil == 1, serotype], "VT", "NVT"), st_index = st_index)] %>%
      glm(st_index == "sero1_final" ~ sttype, data=., binomial)
    x %>% parseGLM("sttype") %>% .[, variable := "Serotype dominant (all participants)"] %>% .[]
  },
  {
    x = lab_data %>%
      .[pneusil_carr == 1 & nvtpneusil_carr == 1] %>%
      melt(measure.vars = c("sero1_final", "sero2_final", "sero3_final", "sero4_final", "sero5_final"),
           variable.name="st_index", value.name="st") %>%
      .[st != ""] %>%
      .[, .(sttype = ifelse(st %in% serotypes_by_vaccine[pneumosil == 1, serotype], "VT", "NVT"), st_index = st_index)] %>%
      glm(st_index == "sero1_final" ~ sttype, data=., binomial)
    x %>% parseGLM("sttype") %>% .[, variable := "Serotype is dominant (in those carrying both VT and NVT)"] %>% .[]
  })

tableSC2[, ci_format := paste0(ci_low," - ", ci_high)]
tableSC2[is.na(ci_low) | is.na(ci_high), ci_format := NA_character_]
tableSC2[type == "continuous", c("option", "variable") := .(variable, NA_character_)]
tableSC2_dominant_carriage = tableSC2[, c("variable", "option", "est", "ci_format", "pval", "N")]

if(make_table) tableSC2_dominant_carriage %>%
  kblOut(booktabs=T, align=c("l","l","l", "l","l"), linesep = "",
         col.names=c("Variable", "", "OR", "95% CI", "p-value", "N"), other_functions = list(
           function(x) kable_styling(x, latex_options = "scale_down")), out_name = "tableSC2_dominant_carriage")
  
#' OR of multiple carriage by age - summarized in text
glm(num_sero_final>1 ~ agegrp, family=stats::binomial(link="logit"),
    data=participant_data_design_ps$variables[pneu_carr_final == 1 & !is.na(num_sero_final)]) %>%
  parseGLM("agegrp") %>% .[]