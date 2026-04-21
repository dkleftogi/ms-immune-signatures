# Assisting personalised treatment decisions in multiple sclerosis using data-driven immunological signatures  

R-based workflows for CyTOF preprocessing, harmonization, and immune signature discovery in multiple sclerosis.

Additionally the repo hosts a pipeline to reproduce the analysis presented in our article titled:

"Intrathecal mesenchymal stem cell therapy in progressive multiple sclerosis: cross-compartment immune profiling in the SMART-MS randomized trial"

## Scope
This repository contains analysis code and documentation only.
No sensitive patient-level data are included.


## Relevant publications

1. Intrathecal mesenchymal stem cell therapy in progressive multiple sclerosis: cross-compartment immune profiling in the SMART-MS randomized trial, Under Review to the Journal of Neuroinflammation, 2025

## Funding

The project is funded by the following grants awarded to Dimitrios Kleftogiannis
1. Project title: “Comparative Immune Profiling of Multiple Sclerosis Treatments: Towards Stratified Therapy”, Funder: Gerda Meyer Nyquist Gulbrandson & Gerd Meyer Nyquist legat (2025)

2. Project title: “Pilot study on high-dimensional analysis for personalized treatment decisions in multiple sclerosis”, Funder: MS-forbundet (2025)

3. Project title: “Assisting personalised treatment decisions in multiple sclerosis using data-driven immunological signatures”, Funder: Helse Vest (2024)

## Description of codes and pipelines

1. Folder R/ contains main utility functions and codes for different parts of the analysis

2. Folder scripts/ contains executable pipelines suitable for interactive analysis in RStudio

3. Folder analysis/ contain command-line wrappers as well as markdown versions of the analysis 

## Data availability
This repository contains the public analysis code used in the study. Publicly released data are provided separately through Zenodo. Participant-level clinical covariates used in selected adjusted models are not included in this repository.

These covariates are stored separately because linking pseudonymised participant identifiers with demographic metadata may increase re-identification risk. Pseudonymised data are generally still considered personal data under data protection guidance.

## Reproducibility
The repository is structured so that approved users can reproduce the full analysis by placing the restricted clinical metadata file in:

clinical_covariates_SMART.csv

If the restricted file is not available, parts of the workflow that require participant-level covariates may not run or may be skipped, depending on the script. Please contact the authors for more details. 

## Ethics/privacy note
The public code repository does not include participant-level age/sex metadata linked to study IDs. This separation was chosen to support reproducibility while reducing the risk of indirect identification from combined datasets. Public release therefore includes code and public data resources only, whereas restricted participant-level covariates remain in the approved secure environment.

## Contact

Comments and reports are welcome, please email: Dimitrios Kleftogiannis 
(dimitrios.kleftogiannis@uib.no)

We are also interested to know about how you have used our source codes, including any improvements that you have implemented.
 
You are free to modify, extend or distribute our source code, as long as our copyright notice remains unchanged and included in its entirety. 

### Access to restricted metadata

Access to participant-level clinical covariates is not provided through the public repository. Researchers seeking to reproduce the full adjusted analyses should contact the study authors or the relevant data custodian and obtain any required approvals before access is granted.

## License

This project is licensed under the MIT License.

Copyright 2026 Department of Clinical Medicine (K1), University of Bergen (UiB) and Neuro-SysMed center for clinical treatment research, Norway

You may only use the source code in this repository in compliance with the license provided in this repository. For more details, please refer to the file named "LICENSE.md".

#### History

21-Jan-2026 : 01 part - CyTOF cleanup with filtering and GMMs

4-Feb-2026  : 02 part - CyTOF harmonisation using cyCombine

20-Feb-2026 : 03 part - SMART-MS omics paper analyses 

21-Apr-2026 : SMART-MS repo updates