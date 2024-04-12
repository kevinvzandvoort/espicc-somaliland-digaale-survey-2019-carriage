#############################################################################################################################
#' Section 1. Matched data on participant id
#############################################################################################################################  
matched_pid = household_data_members_match[, c("pid")] %>% copy %>% .[, household_members_data := 1] %>%
  merge(contact_data_match[, c("pid")] %>% copy %>% .[, contact_data := 1], by="pid", all = TRUE) %>%
  merge(nutrition_data_match[, c("pid")] %>% copy %>% .[, nutrition_data := 1], by="pid", all = TRUE) %>%
  merge(swab_data_match[, c("pid", "swab_id")] %>% copy %>% .[, swab_data := 1], by="pid", all = TRUE) %>%
  merge(lab_data_match[, c("swab_id")] %>% copy %>% .[, lab_data := 1], by="swab_id", all = TRUE) %>%
  .[is.na(household_members_data), household_members_data := 0] %>%
  .[is.na(contact_data), contact_data := 0] %>%
  .[is.na(nutrition_data), nutrition_data := 0] %>%
  .[is.na(swab_data), swab_data := 0] %>%
  .[is.na(lab_data), lab_data := 0]

#############################################################################################################################
#' Section 2. Matched data on household id
#############################################################################################################################  
matched_hid = household_data_match[, c("hid")] %>% copy %>% .[, household_data := 1] %>%
  merge(contact_data_match[, c("hid")] %>% copy %>% .[, contact_data := 1], by="hid", all = TRUE) %>%
  merge(nutrition_data_match[, c("hid")] %>% copy %>% .[, nutrition_data := 1], by="hid", all = TRUE) %>%
  merge(swab_data_match[, c("hid", "swab_id")] %>% copy %>% .[, swab_data := 1], by="hid", all = TRUE) %>%
  merge(lab_data_match[, c("swab_id")] %>% copy %>% .[, lab_data := 1], by="swab_id", all = TRUE) %>%
  .[is.na(household_data), household_data := 0] %>%
  .[is.na(contact_data), contact_data := 0] %>%
  .[is.na(nutrition_data), nutrition_data := 0] %>%
  .[is.na(swab_data), swab_data := 0] %>%
  .[is.na(lab_data), lab_data := 0]

#############################################################################################################################
#' Section 3. Matched data on both participant and household id
#############################################################################################################################  
matched_pid_hid = household_data_members_match[, c("pid", "hid")] %>% copy %>% .[, household_members_data := 1] %>%
  merge(contact_data_match[, c("pid", "hid")] %>% copy %>% .[, contact_data := 1], by=c("pid", "hid"), all = TRUE) %>%
  merge(nutrition_data_match[, c("pid", "hid")] %>% copy %>% .[, nutrition_data := 1], by=c("pid", "hid"), all = TRUE) %>%
  merge(swab_data_match[, c("pid", "hid", "swab_id")] %>% copy %>% .[, swab_data := 1], by=c("pid", "hid"), all = TRUE) %>%
  merge(lab_data_match[, c("swab_id")] %>% copy %>% .[, lab_data := 1], by="swab_id", all = TRUE) %>%
  merge(household_data_match[, c("hid")] %>% copy %>% .[, household_data := 1], by="hid", all=TRUE) %>%
  .[is.na(household_members_data), household_members_data := 0] %>%
  .[is.na(household_data), household_data := 0] %>%
  .[is.na(contact_data), contact_data := 0] %>%
  .[is.na(nutrition_data), nutrition_data := 0] %>%
  .[is.na(swab_data), swab_data := 0] %>%
  .[is.na(lab_data), lab_data := 0]

#############################################################################################################################
#' Section 4. Generate Supplemental Figure A1. Matching records between datasets.
#############################################################################################################################  
#matched_pid_plot = matched_pid %>%
# .[!(household_members_data == 1 & (nutrition_data + lab_data + contact_data) == 0), -"swab_data"] %>%
# upset(order.by="freq")
#matched_hid_plot = matched_hid %>% .[!(household_data == 1 & (nutrition_data + lab_data + contact_data) == 0), -"swab_data"] %>%
# upset(order.by="freq")
matched_pid_hid_plot = matched_pid_hid %>%
  .[!(nutrition_data + lab_data + contact_data) == 0, -"swab_data"] %>%
  .[, c("pid", "hid", "household_data", "household_members_data", "contact_data", "nutrition_data", "lab_data")] %>%
  setNames(c("pid", "hid", "Household", "Household members", "Contact", "Nutrition", "Swab")) %>%
  UpSetR::upset(order.by="freq")

ext="png"
png(file=sprintf("%s/output/%s/figures/%s/figureSA1_matched_datasets.%s", analysis_dir, OUTPUT_DIR, ext, ext),
    units = "in", width = plot_double_col, height = plot_single_col, res=300)
print(matched_pid_hid_plot)
dev.off()

ext="tiff"
tiff(filename=sprintf("%s/output/%s/figures/%s/figureSA1_matched_datasets.%s", analysis_dir, OUTPUT_DIR, ext, ext),
    width = plot_double_col, height = plot_single_col, units = "in", res=300)
print(matched_pid_hid_plot)
dev.off()

ext="eps"
setEPS()
postscript(file=sprintf("%s/output/%s/figures/%s/figureSA1_matched_datasets.%s", analysis_dir, OUTPUT_DIR, ext, ext),
           width = plot_double_col, height = plot_single_col)
print(matched_pid_hid_plot)
dev.off()
