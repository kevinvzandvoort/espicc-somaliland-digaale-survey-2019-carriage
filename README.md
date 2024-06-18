# Pre-vaccination carriage prevalence of Streptococcus pneumoniae serotypes among internally displaced people in Somaliland. Data and analyses.

This repository contains the data and code used for all analyses described in our manuscript:

Van Zandvoort K, Hassan AI, Bobe MO, Pell CL, Ahmed MS, Ortika BD, Ibrahim S, Abdi MI, Karim MA, Eggo RM, Ali SY, Hinds J, Soleman MS, Cummings R, McGowan C, Mulholland K, Hergeeye MA, Satzke C, Checchi F, Flasche S, *Pre-vaccination carriage prevalence of Streptococcus pneumoniae serotypes among internally displaced people in Somaliland*. Preprint available at <https://doi.org/10.1101/2024.02.09.24302568>.

This work is part of a larger study: Evaluating Strategies for Pneumococcal Immunization Campaigns in Crises [(ESPICC)](https://www.elrha.org/project/pneumococcal-vaccination-strategies-for-crisis-affected-populations/).

This analysis is closely linked to our previous [manuscript](https://doi.org/10.1016/j.epidem.2022.100625) and [analysis](https://github.com/kevinvzandvoort/espicc-somaliland-digaale-survey-2019) on [Social contacts and other risk factors for respiratory infections among internally displaced people in Somaliland](https://doi.org/10.1016/j.epidem.2022.100625).

### How to download

You can clone the repository or download the zip from this URL: <https://github.com/kevinvzandvoort/espicc-somaliland-digaale-survey-2019-carriage/archive/refs/heads/main.zip>.

### Questionnaires

This survey was implemented using Open Data Kit. Android tablets were provided by LSHTM Open Data Kit <https://opendatakit.lshtm.ac.uk>.

- Questionnaires were programmed as multiple forms in xlsx, and can be found in the `./questionnaire/xlsx` folder in the [espicc-somaliland-digaale-survey-2019](https://github.com/kevinvzandvoort/espicc-somaliland-digaale-survey-2019/tree/main/questionnaire/xls) repository.
- Only the additional swab related form has been added to this repository, as the `./questionnaire/xlsx/s5_swab.xlsx` file.
- Th file was converted to an xls file to work on an ODK server, and used with *ODK Collect* during fieldwork. Unlike some of the other forms, no manual edits were made to this xls file.

### Included data sets

Only a subset of the data collected using all survey forms has been used for this analysis.
Data has been anonymized, and links between individual household members have been removed.
The anonymized data can be used to replicate all analyses, figures, and tables in the manuscript.
All data is stored in the `./data` folder. A data dictionary is provided in `./data/data_dictionary.xlsx`

The following datasets are included:

- `participant_data.csv`
  - Reported household- and individual-level risk-factors
  - collected with the `s1_household`, `s2_contacts`, `s3_anthropometry`, and `s5_swab` forms.
- `lab_data.csv`
  - Results from microbiological analyses
  - Nasopharyngeal swabs assessed using lytA qPCR and serotyped using microarray
- `serotypes_by_vaccine.csv`
  - Serotypes included in different PCV products
- `other_data/*`
  - Pneumococcal carriage prevalence by age from studies conducted in the Gambia, Kenya, Uganda, and Malawi
  - Serotype-specific invasiveness estimates from Løchen et al (2022; https://doi.org/10.1371/journal.pcbi.1009389)
  - Temperature data from our test (s0_test) and pilot (s1_pilot) shipment.
- `poststratification_data/*`
  - Estimated population size of specific strata, from the main `espicc-somaliland-digaale-survey-2019` repository

### Figures and tables

Code for the analysis can be found in the `./scripts` folder.
The analysis can be replicated by running the `index.R` file (in R), which sources these scripts.

Figures and tables will be created in a newly formed `./output` folder
