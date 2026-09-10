# install.packages("dplyr")
# install.packages("tidyr")

library(dplyr)
library(tidyr)

eDNA_merged_Filter_ID_PA <- read.csv("0_1_presence_absence/eDNA_merged_Filter_ID_PA.csv", header = TRUE, sep = ",", dec = ".")
eDNA_merged_Transect_ID_PA <- read.csv("0_1_presence_absence/eDNA_merged_Transect_ID_PA.csv", header = TRUE, sep = ",", dec = ".", check.names = FALSE)
eDNA_merged_MPA_PA <- read.csv("0_1_presence_absence/eDNA_merged_MPA_PA.csv", header = TRUE, sep = ",", dec = ".", check.names = FALSE)
metadata <- read.csv("../data/3_eDNA_cleaning/METADATA_MERGED_edited.csv", header = TRUE, sep = ",", dec = ".")
combined_database <- read.csv("../data/1_species_database/1_2_combined_database_complete.csv", header = TRUE, sep = ",", dec = ".")

### filter for invasive species only 

invasive_eDNA_merged_Filter_ID_PA <- eDNA_merged_Filter_ID_PA %>%
  left_join(
    combined_database %>%
      dplyr::select(species_name_fishbase, invasive),
    by = "species_name_fishbase"
  ) %>%
  filter(invasive == 1)

invasive_eDNA_merged_Transect_ID_PA <- eDNA_merged_Transect_ID_PA %>%
  left_join(
    combined_database %>%
      dplyr::select(species_name_fishbase, invasive, origin),
    by = "species_name_fishbase"
  ) %>%
  filter(invasive == 1)

invasive_eDNA_merged_MPA_PA <- eDNA_merged_MPA_PA %>%
  left_join(
    combined_database %>%
      dplyr::select(species_name_fishbase, invasive, origin),
    by = "species_name_fishbase"
  ) %>%
  filter(invasive == 1)

### richness

## per Filter_ID

invasive_richness_Filter_ID <- invasive_eDNA_merged_Filter_ID_PA %>%
  pivot_longer(
    cols = matches("^(DAR\\.2023\\.\\d{4}|DAR\\.2024\\.\\d{4}|SPY\\d+)$"),
    names_to = "Filter_ID",
    values_to = "reads"
  ) %>%
  filter(reads > 0) %>%
  group_by(Filter_ID) %>%
  summarise(
    invasive_richness_Filter_ID = n_distinct(species_name_fishbase),
    .groups = "drop"
  )

## per Transect_ID

invasive_richness_Transect_ID <- invasive_eDNA_merged_Transect_ID_PA %>%
  pivot_longer(
    cols = matches("^\\d+"),
    names_to = "Transect_ID",
    values_to = "reads"
  ) %>%
  filter(reads > 0) %>%
  group_by(Transect_ID) %>%
  summarise(
    invasive_richness_Transect_ID = n_distinct(species_name_fishbase),
    .groups = "drop"
  )

## per MPA

invasive_richness_MPA <- invasive_eDNA_merged_MPA_PA %>%
  pivot_longer(
    cols = matches("^\\d+"),   
    names_to = "MPA_ID",
    values_to = "reads"
  ) %>%
  filter(reads > 0) %>%
  group_by(MPA_ID) %>%
  summarise(
    invasive_richness_MPA = n_distinct(species_name_fishbase),
    .groups = "drop"
  )

## per transect and origin

# red sea
invasive_richness_Transect_ID_red_sea <- invasive_eDNA_merged_Transect_ID_PA %>%
  filter(origin == "RS") %>%
  pivot_longer(
    cols = matches("^\\d+"),
    names_to = "Transect_ID",
    values_to = "reads"
  ) %>%
  filter(reads > 0) %>%
  group_by(Transect_ID) %>%
  summarise(
    invasive_richness_Transect_ID_red_sea = n_distinct(species_name_fishbase),
    .groups = "drop"
  )

# atlantic
invasive_richness_Transect_ID_atlantic <- invasive_eDNA_merged_Transect_ID_PA %>%
  filter(origin == "A") %>%
  pivot_longer(
    cols = matches("^\\d+"),
    names_to = "Transect_ID",
    values_to = "reads"
  ) %>%
  filter(reads > 0) %>%
  group_by(Transect_ID) %>%
  summarise(
    invasive_richness_Transect_ID_atlantic = n_distinct(species_name_fishbase),
    .groups = "drop"
  )

## per MPA and origin

# red sea
invasive_richness_MPA_red_sea <- invasive_eDNA_merged_MPA_PA %>%
  filter(origin == "RS") %>%
  pivot_longer(
    cols = matches("^\\d+"),
    names_to = "MPA_ID",
    values_to = "reads"
  ) %>%
  filter(reads > 0) %>%
  group_by(MPA_ID) %>%
  summarise(
    invasive_richness_MPA_red_sea = n_distinct(species_name_fishbase),
    .groups = "drop"
  )

# atlantic
invasive_richness_MPA_atlantic <- invasive_eDNA_merged_MPA_PA %>%
  filter(origin == "A") %>%
  pivot_longer(
    cols = matches("^\\d+"),
    names_to = "MPA_ID",
    values_to = "reads"
  ) %>%
  filter(reads > 0) %>%
  group_by(MPA_ID) %>%
  summarise(
    invasive_richness_MPA_atlantic = n_distinct(species_name_fishbase),
    .groups = "drop"
  )

# dir.create("1_richness", showWarnings = FALSE)

write.csv(invasive_richness_Filter_ID , "1_richness/1_1_invasive_richness_Filter_ID.csv", row.names = FALSE)
write.csv(invasive_richness_Transect_ID , "1_richness/1_1_invasive_richness_Transect_ID.csv", row.names = FALSE)
write.csv(invasive_richness_MPA , "1_richness/1_1_invasive_richness_MPA.csv", row.names = FALSE)

write.csv(invasive_richness_Transect_ID_red_sea , "1_richness/1_1_invasive_richness_Transect_ID_red_sea.csv", row.names = FALSE)
write.csv(invasive_richness_Transect_ID_atlantic , "1_richness/1_1_invasive_richness_Transect_ID_atlantic.csv", row.names = FALSE)

write.csv(invasive_richness_MPA_red_sea , "1_richness/1_1_invasive_richness_MPA_red_sea.csv", row.names = FALSE)
write.csv(invasive_richness_MPA_atlantic , "1_richness/1_1_invasive_richness_MPA_atlantic.csv", row.names = FALSE)
