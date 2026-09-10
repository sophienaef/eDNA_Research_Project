# data/3_eDNA_cleaning

## 3_1_eDNA_intial_cleaning.R
step 1 of eDNA cleaning pipeline

## 3_2_eDNA_taxonomic_cleaning.R
step 2 of eDNA cleaning pipeline

## 3_3_eDNA_region_cleaning.R
step 3 of eDNA cleaning pipeline
produces the final cleaned eDNA data


## inputs

### archive_class_ncbi.csv
information on taxonomy of species
used in 3_1_eDNA_intial_cleaning.R

### habitat_genera_taxa_enriched.csv
information on habitat / water types for species
used in 3_2_eDNA_taxonomic_cleaning.R

## unique_genera_only_species_exists.csv
information on genera that have only one species / one species is plausible
used in 3_2_eDNA_taxonomic_cleaning.R

## eDNA_dataset_ACTNOW_edited.csv
eDNA data from ETH workflow

## eDNA_dataset_Montpellier_edited.csv
eDNA data from SPYGEN workflow

## METADATA_MERGED_edited.csv
collected metadata


## outputs
3_1_eDNA_initial_cleaned directory
3_2_eDNA_taxonomic_cleaned directory
3_2_eDNA_region_cleaned directory

