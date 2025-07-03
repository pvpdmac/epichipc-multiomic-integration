# Multi-omic profiling identifies distinct baseline signatures predicting human neonatal antibody responses to hepatitis B vaccine

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Authors:** Casey P. Shannon\*, Kinga K. Smolen\*, Joann Diray-Arce\*, Olubukola T. Idoko\*, Travis M. Blimkie\*, Anita H.J. van den Biggelaar\*, Abhinav K. Checkervarty, Oludare A. Odumade, Andy Y. An, Tue B. Bennike, Rym Ben-Othman, Jing Chen, Bhavjinder K. Dhillon, Alansana Darboe, Reza Falsafi, Rebecca Ford, Annmarie Hoch, Joe Jude, Geraldine Masiria, Kerry McEnaney, Caitlyn McLoughlin, Sebastiano Montante, Elena Morrocchi, Caitlin Syphurs, Arthur Viode, Nelly Amenyogbe, Asimenia Angelidou, Robert Balshaw, Ryan R. Brinkman, Simon D. van Haren, Amy H. Lee, Arnaud Marchant, David Martino, Al Ozonoff, Paolo Palma, William S. Pomat, Peter C. Richmond, Guzman Sanchez-Schmitz, **EPIC-HIPC Consortium**, Robert E.W. Hancock+, Hanno Steen+, Scott J. Tebbutt+, Tobias R. Kollmann+, Beate Kampmann+, Ofer Levy+

\*These authors contributed equally to the development of this manuscript

+These senior authors contributed equally to the development of this manuscript

**Corresponding Authors:**

-   Casey P. Shannon ([Casey.Shannon\@hli.ubc.ca](mailto:Casey.Shannon@hli.ubc.ca){.email})

-   Ofer Levy ([Ofer.Levy\@childrens.harvard.edu](mailto:Ofer.Levy@childrens.harvard.edu){.email})

**Date:** 2025/05/01

------------------------------------------------------------------------

## Table of Contents

1.  [Overview](#1-overview)
2.  [Manuscript](#2-manuscript)
3.  [Repository Structure](#3-repository-structure)
4.  [System Requirements](#4-system-requirements)
5.  [Installation](#5-installation)
6.  [Data](#6-data)
7.  [Workflow / How to Reproduce](#7-workflow--how-to-reproduce)
    -   [Step 1: Data Preprocessing](#step-1-data-preprocessing)
    -   [Step 2: Further Preprocess Validation Samples](#step-2-further-preprocess-validation-samples)
    -   [Step 3: Prepare Cross-validation](#step-3-prepare-cross-validation)
    -   [Step 4: Main Analysis](#step-4-main-analysis)
    -   [Step 5: Figure Generation](#step-5-figure-generation)
8.  [Outputs](#8-outputs)
9.  [License](#9-license)
10. [Citation](#10-citation)

------------------------------------------------------------------------

## 1. Overview {#id_1-overview}

This repository contains the data, code, and figures associated with the manuscript: **"Multi-omic profiling identifies distinct baseline signatures predicting human neonatal antibody responses to hepatitis B vaccine"**. The primary goal of this research was to identify molecular correlates to predict the humoural response to Hepatitis B vaccine in newborns, and to compare and contrast these correlates to those in adults. This repository provides the necessary resources to reproduce the analyses and figures presented in the paper.

## 2. Manuscript {#id_2-manuscript}

-   **Title:** Multi-omic profiling identifies distinct baseline signatures predicting human neonatal antibody responses to hepatitis B vaccine
-   **Status:** Manuscript currently under review.
-   **Link to manuscript:** TBD

**Highlights:**

-   Multi-omic analysis revealed distinct molecular signatures predictive of HBV vaccine responses in newborns with and without maternal antibodies.
-   Pathways associated with robust HBV responses after multiple doses in newborns corresponded to poor immunogenicity in adults.
-   Age-dependent differences in pathways underlying vaccine responses highlight the unique nature of early-life immunity and underscore the importance of age-specific vaccination strategies.

**Keywords:** multi-omic, systems vaccinology, hepatitis B vaccine, newborns

## 3. Repository Structure {#id_3-repository-structure}

A brief overview of the key directories and their contents:

```         
├── data/
│ ├── raw/                   # Immport data goes here
│ ├── processed/             # Intermediate, processed data files
│ └── external/              # GEO datasets (adult studies)
├── code/
│ ├── 01_data_preprocessing/ # Clean and prepare data
│ ├── 02_cross_validation/   # Run kfoldcv, fit final models
│ ├── 03_figure_generation/  # Generate figures
│ └── utils/                 # Helper functions
├── results/
│ ├── figures/               # Figures
│ ├── tables/                # Tables
│ └── other_outputs/         # Other results, model outputs, etc.
├── requirements.txt         # Environment/dependency specification
├── download_immport_data.md # Instructions to obtain raw data files
├── LICENSE                  # License file for the code
└── README.md                # This file
```

## 4. System Requirements {#id_4-system-requirements}

Large muilti-omic datasets were stored on Amazon Web Services (AWS) S3 for scalable, secure data storage and efficient access. The analysis code was developed and run on an AWS EC2 r5.8xlarge instance, but the code that generates the final figures from intermediate outputs should run on a more modest r5.xlarge instance.

-   **Operating system:** Ubuntu 24.04.1 LTS
-   **Programming language:** R version 4.4.2
-   **Key packages and software:**
    -   `mixOmics` 6.30.0
    -   `fgsea` 1.32.0
    -   `nlme` 3.1-66
    -   `tidyverse` 2.0.0
    -   `quarto` 1.4.555
-   **Hardware requirements:**
    -   Cross-validation: `r5.8xlarge` or equivalent (32 vCPU cores, 256GB RAM)
    -   Figures: `r5.xlarge` or equivalent (4 vCPU cores, 32GB RAM)

## 5. Installation {#id_5-installation}

We deployed the Posit Workbench AMI, which is based on Ubuntu Linux (24.04.1 LTS) and includes software tools such as Python, R, RStudio Server, and Quarto. The instructions below assume a similar starting point (e.g., using one of Docker images provided by [Posit](https://hub.docker.com/u/rstudio).

1.  **Clone the repository:**

    ``` bash
    # bash
    git clone https://github.com/pvpdmac/epichipc/tree/main/Multiomic_integration
    cd Multiomic_integration
    ```

2.  **Set up the computational environment:**

    Install necessary R packages.

    ``` r
    # R
    install.packages('pak')
    pak::pak(readLines('requirements.txt'))
    ```

    Install support for in-line `svg` images (can be conveniently saved from the reports with a right-click).

    ``` bash
    # bash
    quarto install extension coursekata/quarto-inline-svg
    ```

3.  **Download large primary data files:** Refer to `download_immport_data.md`

## 6. Data {#id_6-data}

-   **Source:** 10.3389/fped.2020.00197

-   **Description:** This dataset contains multi-omic data from two neonatal vaccination cohorts:

**Primary cohort (Gambia):** 720 neonates randomized to receive HBV, BCG, or HBV+BCG vaccines at birth (day 0), with catch-up vaccination (BCG+OPV, HBV+OPV, or OPV, respectively) at day 1, 3, or 7, or delayed vaccination (HBV+BCG+OPV) at day 1, 3, or 7.

**Validation cohort (Papua New Guinea):** 101 neonates randomized to the same vaccination groups but only receiving catch-up vaccination at day 7.

**Samples:** blood samples were collected at baseline (day 0) and prior to vaccination on days 1, 3, or 7. Multi-omic profiling included: transcriptomics, epigenetics, flow cytometry, proteomics, and metabolomics.

**Application:** the dataset supported development of machine learning models to predict humoral responce to Hepatitis B vaccine, validated across the PNG cohort, Gambian test set, and *in vitro* tissue constructs.

-   **Raw data:** The data matrices obtained from [Immport](https://immport.org/home "The Immunology Database and Analysis Portal") or [GEO](https://www.ncbi.nlm.nih.gov/geo/ "The Gene Expression Omnibus") should be placed in:

    -   SDY1538:`data/raw/gam/`
    -   SDY2584:`data/raw/png/`
    -   SDY2312:`data/raw/invitro/`

-   **Preprocessing:**

    -   Scripts under `code/01_data_preprocessing/gam/` transform raw data into the analysis-ready format stored in `data/processed/gam/`
    -   Scripts under `code/01_data_preprocessing/png/` transform raw data into the analysis-ready format stored in `data/processed/png/`
    -   `code/01_data_preprocessing/prep_validation/validation_wrangle_gam.R` identifies samples from the Gambia with incomplete omic profiles to serve as an indenpendent test cohort.
    -   `code/01_data_preprocessing/prep_validation/validation_wrangle_png.R` carries out further processing on samples from our validation cohort (assessing feature overlap between gam and png, batch-correction, etc.)

-   **Data dictionary:**

    -   `data/raw/gam/data_disctionary_gam.csv`
    -   `data/raw/png/data_disctionary_png.csv`

## 7. Workflow / How to Reproduce {#id_7-workflow--how-to-reproduce}

To reproduce the analysis and figures, follow these steps in order. All scripts should be run from the root directory of this repository.

### Step 1: Data Preprocessing {#step-1-data-preprocessing}

-   **Purpose:** To clean raw data, merge datasets, perform transformations, and prepare data for analysis. Scripts should be run in the order presented below.

-   **Preprocess GAM data:**

    -   `code/01_data_preprocessing/gam/import_metadata.R`
    -   `code/01_data_preprocessing/gam/import_responses.R`
    -   `code/01_data_preprocessing/gam/import_flow_flowtype.R`
    -   `code/01_data_preprocessing/gam/import_metabolomics.R`
    -   `code/01_data_preprocessing/gam/import_luminex_cytokines.R`
    -   `code/01_data_preprocessing/gam/import_proteomics.R`
    -   `code/01_data_preprocessing/gam/import_transcriptomics.R`
    -   `code/01_data_preprocessing/gam/import_epigenetics.R`
    -   `code/01_data_preprocessing/gam/process_exclusion_list.R`

-   **Preprocess PNG data:**

    -   `code/01_data_preprocessing/png/import_metadata.R`
    -   `code/01_data_preprocessing/png/import_responses.R`
    -   `code/01_data_preprocessing/png/import_luminex_cytokines.R`
    -   `code/01_data_preprocessing/png/import_proteomics.R`
    -   `code/01_data_preprocessing/png/import_transcriptomics.R`
    -   `code/01_data_preprocessing/png/import_epigenetics.R`

-   **Prepare Validation Sets:**

    -   `code/01_data_preprocessing/process_validation/validation_wrangle_gam.R`
    -   `code/01_data_preprocessing/process_validation/validation_wrangle_png.R`

-   **To run:**

    ``` bash
    # bash
    Rscript code/01_data_preprocessing/import_metadata.R
    ...
    ```

-   **Run time:** \<30min

-   **Inputs:** Files from `data/raw/...`

-   **Outputs:** Processed data files in `data/processed/...` (e.g., `data/processed/gam/import_metadata.rds`)

### Step 2: Further Preprocess Validation Samples {#step-2-further-preprocess-validation-samples}

-   **Purpose:** To set aside validation samples, address technical batch-effects between study cohorts.

-   **Prepare Validation Sets:**

    -   `code/01_data_preprocessing/process_validation/validation_wrangle_gam.R`
    -   `code/01_data_preprocessing/process_validation/validation_wrangle_png.R`

-   **To run:**

    ``` bash
    # bash
    Rscript code/01_data_preprocessing/process_validation/validation_wrangle_gam.R
    ...
    ```

-   **Run time:** \<5min

-   **Inputs:** Files from `data/processed/...`

-   **Outputs:** Processed data files in `data/processed/process_validation/...`

### Step 3: Prepare Cross-validation {#step-3-prepare-cross-validation}

-   **Purpose:** To create tibbles to organize data subsets and hyperparameter grid search.

-   **Prepare cross-validation tibble:**

    -   `code/01_data_preprocessing/cross_validation/01_create_data_master_set.R`
    -   `code/01_data_preprocessing/cross_validation/02_create_glmnet_queries_master_set.R`
    -   `code/01_data_preprocessing/cross_validation/03_create_diablo_queries_master_set.R`

-   **To run:**

    ``` bash
    # bash
    Rscript code/01_data_preprocessing/cross_validation/01_create_data_master_set.R
    ...
    ```

-   **Run time:** \<10min

-   **Inputs:** Files from `data/processed/...`

-   **Outputs:** Processed data files in `data/processed/cross_validation/...` (e.g., `data/processed/cross_validation/data_master_set_top20p.rds`)

### Step 4: Main Analysis {#step-4-main-analysis}

-   **Purpose:** Carry out cross-validation, summarize performance, fit final models, etc.

-   **Scripts:**

    -   `code/02_cross_validation/01_batch_run_singleomics.R`
    -   `code/02_cross_validation/02_batch_performance_singleomics.R`
    -   `code/02_cross_validation/03_batch_run_multiomics.R`
    -   `code/02_cross_validation/04_batch_performance_multiomics.R`

-   **To run:**

    ``` bash
    # bash
    Rscript code/02_cross_validation/01_batch_run_singleomics.R
    ...
    ```

-   **Run time:** the cross-validation scripts are likely too slow to replicate (\>20d on the recommended hardware).

-   **Inputs:** Files from `data/processed/cross_validation/...`

-   **Outputs:** Model objects, statistical summaries, tables as serialized R objects in `results/other_output/...`.

### Step 5: Figure Generation {#step-5-figure-generation}

-   **Purpose:** To generate the figures presented in the manuscript.

-   **Scripts:**

    -   `code/03_figure_generation/figure_2_S3.qmd`
    -   `code/03_figure_generation/figure_3_S5_S6.qmd`
    -   `code/03_figure_generation/figure_4.qmd`
    -   `code/03_figure_generation/figure_5_S8.qmd`
    -   `code/03_figure_generation/figure_S4.qmd`
    -   `code/03_figure_generation/figure_S7.qmd`
    -   `code/03_figure_generation/figure_S9.qmd`

-   **To run:**

    ``` bash
    # bash
    quarto render code/03_figure_generation/figure_2_S3.qmd
     ...
    ```

-   **Run time:** \<10min

-   **Inputs:** Files from `data/processed/...` and `results/other_outputs/...` (e.g., saved model objects).

-   **Outputs:** Images (in `.svg` format) embedded in `.html` files, found in `results/figures/...`.

## 8. Outputs {#id_8-outputs}

-   **Figures:** All manuscript figures are located in `results/figures/...`.
    -   `results/figures/figure_2_S3.html`:
        -   Fig. 2: HBV HBsAb responses
        -   Fig. S3: Effect of maternal HBsAb titers
    -   `results/figures/figure_3_S5_S6.html`:
        -   Fig. 3: Multi-omic model performance
        -   Fig. S5: Performance when applied to newborns with discordant HBsAb status at birth
        -   Fig. S6: Comparison of HBV HBsAb responses between study sites
    -   `results/figures/figure_4.html`:
        -   Fig. 4: Pathways implicated by multi-omic model
    -   `results/figures/figure_5_S8.html`:
        -   Fig. 5, S8: Comparison with molecular correlates of HBV response in adults
    -   `results/figures/figure_S4.html`:
        -   Fig. S4: Overview of cross-validation results
    -   `results/figures/figure_S7.html`:
        -   Fig. S7: Effect of vaccination with HBV, BCG, HBV+BCG on molecular correlates *in vitro*
    -   `results/figures/figure_S9.html`:
        -   Fig. S9: Comparison of shared coefficient weights between multi-omic models
-   **Other important outputs:**
    -   `results/other_outputs/trained_model.rds`:
        -   Final models to predict HBV humoral response from a day of birth blood draw
        -   Models are stored in the column `diablo`. These are standard `mixOmics::block.spls` objects.

## 9. License {#id_9-license}

The code in this repository is licensed under the MIT License. See the `LICENSE` file for more details.

Data used in this project may be subject to different licensing or usage terms; please refer to the original data sources for this information.

## 10. Citation {#id_10-citation}

If you use the code or data from this repository in your research, please cite:

-   **The manuscript:** @TODO

-   **This repository:** @TODO

------------------------------------------------------------------------
