# Instructions to download analysis data from ImmPort

1. Register for ImmPort account, if necessary
2. Login to ImmPort
3. Navigate to SDY1538 (GAM) and SDY2584 (PNG) under Shared Data
4. Install IBM Aspera Connect
5. Download the analysis data files needed
6. Save data files to /Downloads folder on local computer to read in for analysis

## Register for ImmPort account, if necessary
In order to download shared data from ImmPort, you will first need to register for an account. To do so, navigate to https://www.immport.org/auth/login 

<img width="1387" alt="Image" src="https://github.com/user-attachments/assets/5f33e2ee-bc6a-482e-a868-4967e968c86a" />

## Login to ImmPort
Enter user credentials to login to ImmPort account.

<img width="1387" alt="Image" src="https://github.com/user-attachments/assets/6d9230b6-8bc1-4e28-a6fd-8d2a03ea5fdc" />

## Navigate to SDY1538 (GAM) and SDY2584 (PNG)

**Gambia Main (GAM)**

- To access shared data for the Gambia Main (GAM) cohort, please navigate to [SDY1538](https://www.immport.org/shared/study/SDY1538) 

- For GAM epigenetics data files, please refer to GEO accession [GSE272800](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272800)

**Papua New Guinea (PNG)**

- To access shared data for the Papua New Guinea (PNG) cohort, please navigate to [SDY2584](https://www.immport.org/shared/study/SDY2584) 

- For PNG epigenetics data files, please refer to GEO accession [GSE270978](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270978)

## Install IBM Aspera Connect
The IBM Aspera Connect browser extension and desktop software must be installed before downloading data from ImmPort.

<img width="1405" alt="Image" src="https://github.com/user-attachments/assets/01fc4c9d-9121-4f08-bf34-71dee936a8a1" />

## Download the analysis data files
You will need to download the following list of data files from ImmPort. Navigate to the 'Study Files' tab under the corresponding SDY, select the files listed, and select the Download button.

SDY1538 (GAM) files:
- GAMMAIN_SINGLEOMICS_ADA.csv
- GAMMAIN_SINGLEOMICS_BCELLFREQS.csv
- GAMMAIN_SINGLEOMICS_CLINICAL.csv
- GAMMAIN_SINGLEOMICS_CYT.csv
- GAMMAIN_SINGLEOMICS_MET_COUNTS_WITHXENO.csv
- GAMMAIN_SINGLEOMICS_MYELOIDFREQS.csv
- GAMMAIN_SINGLEOMICS_PRT.csv
- GAMMAIN_SINGLEOMICS_RNA_COUNTS.csv
- GAMMAIN_SINGLEOMICS_TITERS.csv
- GAM_Main_2020-05-01_EssentialClinicalData_Clean.csv
- GAM_Main_Tier1_Clean.csv
- GAM_RandomizationGroupDefinition.csv
- exclusion_list.csv

SDY2584 (PNG) files:
- PNGVAL_SINGLEOMICS_BCELLS_FREQS.txt
- PNGVAL_SINGLEOMICS_MYELOID_FREQS.txt
- PNGVAL_SINGLEOMICS_CLINICAL.csv
- PNGVAL_SINGLEOMICS_TITERS.csv
- PNGVAL_SINGLEOMICS_CYT_invivo_updated_15Dec2023.csv
- PNGVAL_SINGLEOMICS_RNA_COUNTS_V1V2.csv
- PNGVAL_SINGLEOMICS_PRT_COUNTS_V1V2.csv
- PNGVAL_SINGLEOMICS_MET_COUNTS_V1V2.csv
- PNGVAL_SINGLEOMICS_ADA_invivo.csv

<img width="1118" alt="Image" src="https://github.com/user-attachments/assets/34a938c7-1804-4ab5-b3e3-14d446d37e9d" />

## Save files to /Downloads folder on local computer to read in for analysis
Aspera Connect will download the files to your local /Downloads folder. You will need to manually update all file directories within the analysis code to these local directories. 

<img width="491" alt="Image" src="https://github.com/user-attachments/assets/c836047b-b6b6-491d-90c9-14b7da83d379" />
