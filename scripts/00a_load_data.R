#' participant data
#' - includes data collected during the household and participant survey
participant_data = fread("./data/participant_data.csv")
#participant_data = participant_data[, -c("animal_proximity", "animal_touch", "animal_bitten", 
#                      "animal_kill", "animal_cook")]
#participant_data = participant_data[, -c("travel_furthest_name", "travel_furthest_name_other")]
#participant_data = participant_data[, -c("nutrition_issue", "nutrition_issue_equipment")]
#fwrite(participant_data, "./data/participant_data.csv")

#' lab data
#' - additional data for each collected swab
lab_data = fread("./data/lab_data.csv")
#lab_data = lab_data[, -"hid"]
#fwrite(lab_data, "./data/lab_data.csv")

#' serotypes included in each PCV product
serotypes_by_vaccine = data.table::fread("./data/serotypes_by_vaccine.csv")

#' Nb any data that may link individuals to the same household has been removed to ensure participant data
#'  remains non-identifiable.