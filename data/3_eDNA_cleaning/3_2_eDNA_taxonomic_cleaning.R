# install.packages("dplyr")
# install.packages("rfishbase")

library(dplyr)
library(rfishbase)

eDNA_merged <- read.csv("3_1_eDNA_initial_cleaned/3_1_eDNA_merged.csv", header = TRUE, sep = ",", dec = ".")
habitat_genera_taxa_enriched <- read.csv("habitat_genera_taxa_enriched.csv", header = TRUE, sep = ",", dec = ".")
unique_genera_only_species_exists <- read.csv("unique_genera_only_species_exists.csv", header = TRUE, sep = ",", dec = ".")

### delete rows with 0 reads

eDNA_merged <- eDNA_merged %>%
  filter(!if_all(12:ncol(eDNA_merged), ~ .x == 0))

### taxonomic cleaning

# filter out any taxonomy that isn't defined to the family level

allowed_ranks <- c("subspecies", "species", "genus", "family", "subfamily", "tribe")
eDNA_merged <- eDNA_merged %>%
  filter(rank %in% allowed_ranks)

# remove families where no lower classification was found

filtered_families <- eDNA_merged %>%
  filter(
    rank == "family",
    is.na(genus_name),
    is.na(species_name)
  )

unique_filtered_families <- filtered_families %>%
  distinct(family_name) %>%
  pull(family_name)

families_with_lower <- eDNA_merged %>%
  filter(rank %in% c("genus", "species", "subspecies")) %>%
  distinct(family_name) %>%
  pull(family_name)

eDNA_merged <- eDNA_merged %>%
  filter(
    !(
      rank == "family" &
        family_name %in% unique_filtered_families &
        !(family_name %in% families_with_lower)
    )
  )

### demote and clean certain ranks 

# Notoscopelus elongatus kroyeri -subspecies

eDNA_merged <- eDNA_merged %>%
  mutate(
    scientific_name = ifelse(
      species_name == "Notoscopelus elongatus" & rank == "subspecies",
      "Notoscopelus elongatus",
      scientific_name
    ),
    rank = ifelse(
      species_name == "Notoscopelus elongatus" & rank == "subspecies",
      "species",
      rank
    )
  )

# Scombridae - Thunnini tribe

# Xenocyprididae - Xenocypridinae	subfamily
# Cobitidae - Cobitinae	subfamily
# Salmonidae - Salmoninae subfamily
# Leuciscidae _ Laviniinae and Leuciscinae subfamily

targets_subfamily <- c(
  "Xenocypridinae",
  "Cobitinae",
  "Laviniinae",
  "Leuciscinae",
  "Salmoninae"
)
targets_tribe <- c(
  "Thunnini"
)

eDNA_merged <- eDNA_merged %>%
  mutate(
    change_to_family = 
      (scientific_name %in% targets_tribe & rank == "tribe") |
      (scientific_name %in% targets_subfamily & rank == "subfamily"),
    
    scientific_name = ifelse(change_to_family, family_name, scientific_name),
    rank = ifelse(change_to_family, "family", rank)
  ) %>%
  select(-change_to_family)

### remove freshwater 

## species

# add fishbase valid species names

valid_species_lookup <- eDNA_merged %>%
  filter(!is.na(species_name)) %>%
  distinct(species_name)
validated <- validate_names(valid_species_lookup$species_name)
valid_species_lookup$species_name_fishbase <- validated

# adjust some manually

valid_species_lookup <- valid_species_lookup %>%
  mutate(
    species_name_fishbase = case_when(
      species_name == "Pegusa bleekeri" ~ "Ebosia bleekeri",
      species_name == "Repomucenus filamentosus" ~ "Callionymus filamentosus",
      TRUE ~ species_name_fishbase
    )
  )

# merge 

eDNA_merged <- eDNA_merged %>%
  left_join(
    valid_species_lookup %>% select(species_name, species_name_fishbase),
    by = "species_name"
  ) %>%
  relocate(species_name_fishbase, .after = 9)

# add water type

fish_info <- species(eDNA_merged$species_name_fishbase)

fish_info <- fish_info %>%
  select(
    Species,
    Fresh,
    Brack,
    Saltwater
  ) %>%
  rename(species_name_fishbase = Species)

eDNA_merged <-
  eDNA_merged %>%
  left_join(fish_info, by = "species_name_fishbase") %>%
  relocate(names(fish_info)[names(fish_info) != "species_name_fishbase"], .after = 11)

eDNA_merged <- eDNA_merged %>%
  left_join(
    habitat_genera_taxa_enriched %>%
      select(taxon_name, Fresh_taxa = Fresh, Brack_taxa = Brack, Saltwater_taxa = Saltwater),
    by = c("scientific_name" = "taxon_name")
  ) %>%
  mutate(
    Fresh = ifelse(is.na(Fresh), Fresh_taxa, Fresh),
    Brack = ifelse(is.na(Brack), Brack_taxa, Brack),
    Saltwater = ifelse(is.na(Saltwater), Saltwater_taxa, Saltwater)
  ) %>%
  select(-Fresh_taxa, -Brack_taxa, -Saltwater_taxa)

## remove where saltwater == 0

eDNA_merged <- eDNA_merged %>%
  filter(Saltwater != 0) %>%
  select(-Fresh, -Brack, -Saltwater)

### complete 0.98 threshold

## remove all non-genera where best_identity is < 0.98 

eDNA_merged <- eDNA_merged %>%
  filter(!(rank != "genus" & best_identity < 0.98))

## fill in species if only one exists in genus

eDNA_merged <- eDNA_merged %>%
  left_join(
    unique_genera_only_species_exists %>%
      select(genus_name, species_name),
    by = "genus_name",
    suffix = c("", "_lookup")
  ) %>%
  mutate(
    species_name = if_else(
      rank == "genus" & !is.na(species_name_lookup),
      species_name_lookup,
      species_name
    ),
    scientific_name = if_else(
      rank == "genus" & !is.na(species_name_lookup),
      species_name_lookup,
      scientific_name
    ),
    rank = if_else(
      rank == "genus" & !is.na(species_name_lookup),
      "species",
      rank
    )
  ) %>%
  select(-species_name_lookup)

# add fishbase valid species names for new ones

valid_species_lookup <- eDNA_merged %>%
  filter(!is.na(species_name), is.na(species_name_fishbase)) %>%
  distinct(species_name)
validated <- validate_names(valid_species_lookup$species_name)
valid_species_lookup$species_name_fishbase <- validated

eDNA_merged <- eDNA_merged %>%
  left_join(
    valid_species_lookup,
    by = "species_name",
    suffix = c("", "_new")
  ) %>%
  mutate(
    species_name_fishbase = coalesce(species_name_fishbase, species_name_fishbase_new)
  ) %>%
  select(-species_name_fishbase_new)

### change definitely wrong species

eDNA_merged <- eDNA_merged %>%
  mutate(
    species_name = if_else(scientific_name == "Labrus mixtus",
                           "Labrus bimaculatus",
                           species_name),
    species_name_fishbase = if_else(scientific_name == "Labrus mixtus",
                                    "Cichlasoma bimaculatum",
                                    species_name_fishbase),
    scientific_name = if_else(scientific_name == "Labrus mixtus",
                              "Labrus bimaculatus",
                              scientific_name),
  )

# add fishbase genus

eDNA_merged <- eDNA_merged %>%
  mutate(
    genus_name_fishbase = word(species_name_fishbase, 1),
    genus_name_fishbase = coalesce(genus_name_fishbase, genus_name)
  ) %>%
  relocate(genus_name_fishbase, .after = 8)

### write files

# dir.create("3_2_eDNA_taxonomic_cleaned", showWarnings = FALSE)

write.csv(eDNA_merged, "3_2_eDNA_taxonomic_cleaned/3_2_eDNA_merged.csv", row.names = FALSE)
