# install.packages("dplyr")
# install.packages("tidyr")

library(dplyr)
library(tidyr)

eDNA_merged <- read.csv("3_2_eDNA_taxonomic_cleaned/3_2_eDNA_merged.csv", header = TRUE, sep = ",", dec = ".")
metadata <- read.csv("METADATA_MERGED_edited.csv", header = TRUE, sep = ",", dec = ".")

### remove technical columns and reorder

eDNA_merged <- eDNA_merged[, -c(1:6)]
eDNA_merged <- eDNA_merged %>%
  select(6, 7, 4, 5, 2, 3, 1, everything())

### pool reads per taxa

eDNA_merged <- eDNA_merged %>%
  group_by(
    scientific_name,
    rank,
    species_name,
    species_name_fishbase,
    genus_name,
    genus_name_fishbase,
    family_name
  ) %>%
  summarise(
    across(everything(), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

### invasive genera broadening

eDNA_merged <- eDNA_merged %>%
  mutate(
    species_name = case_when(
      genus_name %in% c("Hemiramphus", "Torquigener") ~
        case_when(
          genus_name == "Hemiramphus" ~ "Hemiramphus far",
          genus_name == "Torquigener" ~ "Torquigener flavimaculosus"
        ),
      TRUE ~ species_name
    ),
    species_name_fishbase = case_when(
      genus_name %in% c("Hemiramphus", "Torquigener") ~
        case_when(
          genus_name == "Hemiramphus" ~ "Hemiramphus far",
          genus_name == "Torquigener" ~ "Torquigener flavimaculosus"
        ),
      TRUE ~ species_name_fishbase
    ),
    rank = if_else(
      genus_name %in% c("Hemiramphus", "Torquigener"),
      "species",
      rank
    ),
    scientific_name = if_else(
      genus_name %in% c("Hemiramphus", "Torquigener"),
      species_name,
      scientific_name
    )
  )

#####

### sort eDNA reads per Filter_ID

eDNA_merged_Filter_ID <- eDNA_merged

id_lookup <- setNames(metadata$Filter_ID, metadata$PCR_ID)

# rename columns
colnames(eDNA_merged_Filter_ID)[7:ncol(eDNA_merged_Filter_ID)] <- sapply(
  colnames(eDNA_merged_Filter_ID)[7:ncol(eDNA_merged_Filter_ID)],
  function(x) {
    
    # keep SPY columns unchanged
    if (startsWith(x, "SPY")) {
      return(x)
    }
    
    # replace PCR_ID with Filter_ID
    if (x %in% names(id_lookup)) {
      return(id_lookup[[x]])
    } else {
      return(x)
    }
  }
)

# sort 
eDNA_merged_Filter_ID <- eDNA_merged_Filter_ID %>%
  select(1:7, sort(names(.)[8:ncol(.)]))

# dir.create("3_3_eDNA_region_cleaned", showWarnings = FALSE)

write.csv(eDNA_merged_Filter_ID, "3_3_eDNA_region_cleaned/3_3_eDNA_merged_Filter_ID.csv", row.names = FALSE)


### pool reads 

read_cols <- intersect(colnames(eDNA_merged), metadata$PCR_ID)

## per Transect

eDNA_merged_Transect_ID <- eDNA_merged %>%
  pivot_longer(
    cols = all_of(read_cols),
    names_to = "PCR_ID",
    values_to = "reads"
  ) %>%
  left_join(
    metadata %>% select(PCR_ID, Transect_ID),
    by = "PCR_ID"
  ) %>%
  group_by(across(-c(PCR_ID, reads, Transect_ID)), Transect_ID) %>%
  summarise(
    reads = sum(reads, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Transect_ID,
    values_from = reads,
    values_fill = 0
  )

# dir.create("3_3_eDNA_region_cleaned", showWarnings = FALSE)

write.csv(eDNA_merged_Transect_ID, "3_3_eDNA_region_cleaned/3_3_eDNA_merged_Transect_ID.csv", row.names = FALSE)

## per MPA

eDNA_merged_MPA <- eDNA_merged %>%
  pivot_longer(
    cols = all_of(read_cols),
    names_to = "PCR_ID",
    values_to = "reads"
  ) %>%
  left_join(
    metadata %>% select(PCR_ID, MPA_Identi),
    by = "PCR_ID"
  ) %>%
  group_by(across(-c(PCR_ID, reads, MPA_Identi)), MPA_Identi) %>%
  summarise(
    reads = sum(reads, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = MPA_Identi,
    values_from = reads,
    values_fill = 0
  )

# dir.create("3_3_eDNA_region_cleaned", showWarnings = FALSE)

write.csv(eDNA_merged_MPA, "3_3_eDNA_region_cleaned/3_3_eDNA_merged_MPA.csv", row.names = FALSE)


## by ecoregion

metadata$Ecoregion <- gsub(" ", "_", metadata$Ecoregion)

eDNA_merged_Ecoregion <- eDNA_merged %>%
  pivot_longer(
    cols = all_of(read_cols),
    names_to = "PCR_ID",
    values_to = "reads"
  ) %>%
  left_join(
    metadata %>% select(PCR_ID, Ecoregion),
    by = "PCR_ID"
  ) %>%
  group_by(across(-c(PCR_ID, reads, Ecoregion)), Ecoregion) %>%
  summarise(
    reads = sum(reads, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Ecoregion,
    values_from = reads,
    values_fill = 0
  )

write.csv(eDNA_merged_Ecoregion, "3_3_eDNA_region_cleaned/3_3_eDNA_merged_Ecoregion.csv", row.names = FALSE)
