# install.packages("dplyr")
# install.packages("tidyr")

library(dplyr)
library(tidyr)

eDNA_merged_Filter_ID <- read.csv("../data/3_eDNA_cleaning/3_3_eDNA_region_cleaned/3_3_eDNA_merged_Filter_ID.csv", header = TRUE, sep = ",", dec = ".")
eDNA_merged_Transect_ID <- read.csv("../data/3_eDNA_cleaning/3_3_eDNA_region_cleaned/3_3_eDNA_merged_Transect_ID.csv", header = TRUE, sep = ",", dec = ".", check.names = FALSE)
eDNA_merged_MPA <- read.csv("../data/3_eDNA_cleaning/3_3_eDNA_region_cleaned/3_3_eDNA_merged_MPA.csv", header = TRUE, sep = ",", dec = ".", check.names = FALSE)
eDNA_merged_Ecoregion <- read.csv("../data/3_eDNA_cleaning/3_3_eDNA_region_cleaned/3_3_eDNA_merged_Ecoregion.csv", header = TRUE, sep = ",", dec = ".")

# ### calc. total number of species and genera detected
# unique_species <- eDNA_merged_Filter_ID %>%
#   summarise(unique_fishbase_species = n_distinct(species_name_fishbase, na.rm = TRUE))
# 
# unique_genera <- eDNA_merged_Filter_ID %>%
#   summarise(unique_fishbase_species = n_distinct(genus_name_fishbase, na.rm = TRUE))

### convert to presence / absence

eDNA_merged_Filter_ID_PA <- eDNA_merged_Filter_ID %>%
  mutate(
    across(
      8:ncol(.),
      ~ ifelse(. > 0, 1, 0)
    )
  )

eDNA_merged_Transect_ID_PA <- eDNA_merged_Transect_ID %>%
  mutate(
    across(
      8:ncol(.),
      ~ ifelse(. > 0, 1, 0)
    )
  )
eDNA_merged_MPA_PA <- eDNA_merged_MPA %>%
  mutate(
    across(
      8:ncol(.),
      ~ ifelse(. > 0, 1, 0)
    )
  )


eDNA_merged_Ecoregion_PA <- eDNA_merged_Ecoregion %>%
  mutate(
    across(
      8:ncol(.),
      ~ ifelse(. > 0, 1, 0)
    )
  )

### write files

# dir.create("0_1_presence_absence", showWarnings = FALSE)

write.csv(eDNA_merged_Filter_ID_PA , "0_1_presence_absence/eDNA_merged_Filter_ID_PA.csv", row.names = FALSE)
write.csv(eDNA_merged_Transect_ID_PA, "0_1_presence_absence/eDNA_merged_Transect_ID_PA.csv", row.names = FALSE)
write.csv(eDNA_merged_MPA_PA, "0_1_presence_absence/eDNA_merged_MPA_PA.csv", row.names = FALSE)
write.csv(eDNA_merged_Ecoregion_PA , "0_1_presence_absence/eDNA_merged_Ecoregion_PA.csv", row.names = FALSE)
