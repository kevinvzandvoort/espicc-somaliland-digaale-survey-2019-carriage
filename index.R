#' Install and load dependencies
#' - make sure you have pacman installed: install.packages("pacman")
#' - or install and load manually
pacman::p_load(data.table, survey, ggplot2, patchwork, magrittr, kableExtra, stringr, lmtest, zscorer, mitools,
               binom, socialmixr, UpSetR, tiff)

#' Set the analysis directory (usually the root directory of repo)
analysis_dir = getwd()

#' Should tables be created (set to FALSE if correct software is not installed)
make_table = TRUE

#' What variables should be used to generate the post-stratification weights?
#' - see script 0d_poststratification.R for more details
POSTSTRATIFICATION_STRATA = c("age_householdsize", "age_sex", "age", "none")[2]

#' Set global options, create output folders, load helper functions, load data and process data
#' - also creates Survey objects implementing finite population corrections and post-stratification
source(sprintf("%s/scripts/00_setup.R", analysis_dir))

#' Sampling flowchart
#' - generates Figure 1: Flowchart of sampling procedure.
source(sprintf("%s/scripts/01_sampling_flowchart.R", analysis_dir))

#' Descriptive statistics sample and NP swabs
#' - generates Table 1: Sample characteristics and carriage prevalence.
#' - the final section of Table 1, invasive disease likely caused by VT, is generated in the next script
source(sprintf("%s/scripts/02_sample_characteristics.R", analysis_dir))
  
#' Projected invasiveness analysis
#' - generates the final section of Table 1: Invasive disease likely caused by vaccine serotypes.
#' - generates Supplemental Figure D1: Age and serotype specific invasiveness values.
#' - generates Supplemental Figure D2: Estimated proportion of IPD cases caused by serotypes covered by PCVs.
#' - generates Supplemental Table D1: Proportion of current IPD covered by PCVs
#'   Nb. generating the bootstrap samples is not very efficient
source(sprintf("%s/scripts/03_invasiveness_estimates.R", analysis_dir))

#' Serotype distribution
#' - generates Figure 2: Pneumococcal serotype distribution.
#' - generates Supplemental Figure C1. Pneumococcal serotype distribution.
#' - generates Supplemental Figure C2. Pneumococcal serotype distribution (by age).
source(sprintf("%s/scripts/04_serotype_distribution.R", analysis_dir))

#' Carriage prevalence by age
#' - generates Figure 3: Prevalence and serotype distribution by age.
#' - generates Supplemental Figure C3. Prevalence and serotype distribution by age.
source(sprintf("%s/scripts/05_carriage_prevalence.R", analysis_dir))

#' Age-specific exposure to pneumococci
#' - generates Figure 4: The contribution of different age groups towards the age specific exposure to pneumococcus.
#' - generates Supplemental Table C6. The contribution of different age groups towards the age specific exposure to pneumococcus.
#'   Nb. generating the bootstrap samples is not very efficient
source(sprintf("%s/scripts/06_age_specific_exposure.R", analysis_dir))

#' Regression analyses
#' - generates Table 2. Association between risk factors and pneumococcal carriage.
#' - generates Supplemental Table C3. Association between risk factors and pneumococcal density.
source(sprintf("%s/scripts/07_regression_analyses.R", analysis_dir))

#' Other microbiological results
#' - generates Supplemental Table C1. Other microbiological results.
#' - generates Supplemental Table C2. Association between serotype and dominant carriage.
source(sprintf("%s/scripts/08_other_microbiological_results.R", analysis_dir))

#' Sensitivity of population prevalence to post-stratification weights
postratification_sensitivity_data = list()
for(POSTSTRATIFICATION_STRATA in c("age_householdsize", "age_sex", "age", "none")){
  source(sprintf("%s/scripts/00_setup.R", analysis_dir))
  source(sprintf("%s/scripts/09a_poststratification_sensitivity_prevalence.R", analysis_dir))
}
#' - generates Supplemental Table C5. Pneumococcal prevalence estimates by different post-stratification weights
#' - generates Supplemental Figure C4. Prevalence and serotype distribution by age using different weights.
source(sprintf("%s/scripts/09b_poststratification_sensitivity_result.R", analysis_dir))

#' Prevalence compared to other settings
#' - generates Supplemental Figure C5. Prevalence by age compared to different settings.
source(sprintf("%s/scripts/10_other_settings.R", analysis_dir))

#' Matched datasets
#' Nb. To ensure our data remains non-identifiable, any information that links participants to the same household
#'  has been removed from this repository. However, the code used to create the UpSet plot is provided in script 11.
#' - generates Supplemental Figure A1. Matching records between datasets.
#source(sprintf("%s/scripts/11_matched_datasets.R", analysis_dir))

#' Temperature during shipments
#' - generates Supplemental Table B1. Time until temperature exceedance.
#' - generates Supplemental Figure B1. Temperature of test shipment over time.
#' - generates Supplemental Figure B2. Temperature of pilot shipment over time.
source(sprintf("%s/scripts/12_shipment_temperature.R", analysis_dir)) #TODO

#' Prevalence between first (pilot) and second shipment
#' - generates Supplemental Table B2. Carriage prevalence in pilot and second shipment.
source(sprintf("%s/scripts/13_prevalence_shipments.R", analysis_dir))
