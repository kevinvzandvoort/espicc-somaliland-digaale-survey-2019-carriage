#' Values in this script are hardcoded, but equivalent to those in the
#'  https://github.com/kevinvzandvoort/espicc-somaliland-digaale-survey-2019 repository.
#' Download the kevinvzandvoort/espicc-somaliland-digaale-survey-2019 repository to see
#'  calculations for these values.

#' Number of shelters in Digaale
n_shelters = 894

#' Number of shelters where a household member was present / not present when conducting the survey
n_shelters_present = 489
n_shelters_notpresent = 405

#' Number of shelters where a household member consented
n_shelters_consented = 464

#' We conducted a simple random sample of 96 shelters where no individual was present on any visit, and collected
#'  information about 73 shelters
n_shelters_notpresent_sampled = 96
n_shelters_notpresent_sampled_included = 73

#' We use the lower bound of the estimated number of shelters that are vacant to calculate the most conservative
#'  estimate to be used in the fpc
n_shelters_empty = 179

#' We estimate the total number of non-vacant shelters in Digaale
n_shelters_for_fpc = 715

#' We estimate how many of all households living in Digaale were present at the time of the survey, and use the inverse of this number to calculate the fpc_inflate_N_factor
#'  - This estimate will be used to inflate the population size as observed in our study, as an estimate of the total population size living in Digaale
fpc_inflate_N_factor = n_shelters_for_fpc/n_shelters_present