#' participant data
#' - includes data collected during the household and participant survey
participant_data = fread("./data/participant_data.csv")

#' lab data
#' - additional data for each collected swab
lab_data = fread("./data/lab_data.csv")

#' serotypes included in each PCV product
serotypes_by_vaccine = data.table::fread("./data/serotypes_by_vaccine.csv")

#' Nb any data that may link individuals to the same household has been removed to ensure participant data
#'  remains non-identifiable.