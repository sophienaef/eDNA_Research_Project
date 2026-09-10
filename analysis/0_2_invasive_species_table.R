# install.packages("dplyr")
# install.packages("tidyr")

library(dplyr)
library(tidyr)

eDNA_merged_Ecoregion_PA <- read.csv("0_1_presence_absence/eDNA_merged_Ecoregion_PA.csv", header = TRUE, sep = ",", dec = ".")
eDNA_merged_MPA_PA <- read.csv("0_1_presence_absence/eDNA_merged_MPA_PA.csv", header = TRUE, sep = ",", dec = ".", check.names = FALSE)
metadata <- read.csv("../data/3_eDNA_cleaning/METADATA_MERGED_edited.csv", header = TRUE, sep = ",", dec = ".")
combined_database <- read.csv("../data/1_species_database/1_2_combined_database_complete.csv", header = TRUE, sep = ",", dec = ".")

### filter for invasive species only 

invasive_eDNA_merged_MPA_PA <- eDNA_merged_MPA_PA %>%
  left_join(
    combined_database %>%
      dplyr::select(species_name_fishbase, invasive, origin, IUCN_Global),
    by = "species_name_fishbase"
  ) %>%
  filter(invasive == 1)

### summarise

mpa_cols <- names(invasive_eDNA_merged_MPA_PA)[
  grepl("^[0-9]+$", names(invasive_eDNA_merged_MPA_PA))
]


invasive_summary <- invasive_eDNA_merged_MPA_PA %>%
  rowwise() %>%
  mutate(
    number_MPAs = sum(c_across(all_of(mpa_cols)) == 1, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(
    species_name_fishbase,
    family_name,
    origin,
    IUCN_Global,
    number_MPAs
  )%>%
  arrange(species_name_fishbase) 

###

# dir.create("0_2_invasive_species_table", showWarnings = FALSE)

write.csv(invasive_summary , "0_2_invasive_species_table/0_2_invasive_species_table.csv", row.names = FALSE,  quote = FALSE)
