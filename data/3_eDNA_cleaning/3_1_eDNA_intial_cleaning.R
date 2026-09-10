# install.packages("tidyverse")
# install.packages("dplyr")
# install.packages("purrr")

library(tidyverse)
library(dplyr)
library(purrr)

eDNA_Actnow <- read.csv("eDNA_dataset_ACTNOW_edited.csv", header = TRUE, sep = ",", dec = ".")
eDNA_Spygen <- read.csv("eDNA_dataset_Montpellier_edited.csv", header = TRUE, sep = ",", dec = ".")
metadata <- read.csv("METADATA_MERGED_edited.csv", header = TRUE, sep = ",", dec = ".")
archive_class_ncbi <- read.csv("archive_class_ncbi.csv", header = TRUE, sep = ",")

### filter best_identity and assignment

# threshold 0.98 later reinforced

eDNA_Actnow <- eDNA_Actnow[eDNA_Actnow$best_identity >= 0.96, ]
eDNA_Spygen <- eDNA_Spygen[eDNA_Spygen$best_identity >= 0.96, ]

# demote to genus if < 1

condition_actnow <- eDNA_Actnow$best_identity < 1 & !is.na(eDNA_Actnow$species_name)
if (any(condition_actnow)) {
  eDNA_Actnow$species_name[condition_actnow] <- NA
  eDNA_Actnow$rank[condition_actnow] <- "genus"
  eDNA_Actnow$scientific_name[condition_actnow] <- eDNA_Actnow$genus_name[condition_actnow]
}

condition_spygen <- eDNA_Spygen$best_identity < 1 & !is.na(eDNA_Spygen$species_name)
if (any(condition_spygen)) {
  eDNA_Spygen$species_name[condition_spygen]  <- NA
  eDNA_Spygen$rank[condition_spygen] <- "genus"
  eDNA_Spygen$scientific_name[condition_spygen] <- eDNA_Spygen$genus_name[condition_spygen]
}

### add class_name and filter for only Actinopteri and Chondrichthyes

# Actnow
eDNA_Actnow <- eDNA_Actnow %>%
  left_join(
    archive_class_ncbi %>% select(scientific_name, class_name),
    by = "scientific_name"
  ) %>%
  relocate(class_name, .after = 3)

target_taxa <- c(
  "Hemiramphus",
  "Torquigener",
  "Cyprinus",
  "Upeneus",
  "Knipowitschia"
)
eDNA_Actnow <- eDNA_Actnow %>%
  mutate(
    class_name = case_when(
      scientific_name %in% target_taxa & is.na(class_name) ~ "Actinopteri",
      TRUE ~ class_name
    )
  )

eDNA_Actnow <- eDNA_Actnow %>%
  filter(class_name %in% c("Actinopteri", "Chondrichthyes"))

# Spygen
eDNA_Spygen <- eDNA_Spygen %>%
  left_join(
    archive_class_ncbi %>% select(scientific_name, class_name),
    by = "scientific_name"
  ) %>%
  relocate(class_name, .after = 10)

eDNA_Spygen <- eDNA_Spygen %>%
  mutate(
    class_name = case_when(
      order_name %in% c(
        "Anguilliformes",
        "Carangiformes",
        "Characiformes",
        "Cypriniformes",
        "Gobiiformes",
        "Ophidiiformes",
        "Perciformes",
        "Pleuronectiformes",
        "Salmoniformes",
        "Scombriformes",
        "Siluriformes",
        "Stomiiformes"
      ) ~ "Actinopteri",
      
      order_name == "Lamniformes" ~ "Chondrichthyes",
      
      TRUE ~ class_name
    )
  )

eDNA_Spygen <- eDNA_Spygen %>%
  filter(class_name %in% c("Actinopteri", "Chondrichthyes"))

### blanks (Actnow only)

# sort PCR_IDs alphabetically in eDNA_Actnow
eDNA_Actnow <- cbind(
  eDNA_Actnow[, 1:12],
  eDNA_Actnow[, sort(names(eDNA_Actnow)[13:ncol(eDNA_Actnow)])]
)

# PCR_ID named wrong

eDNA_Actnow <- eDNA_Actnow %>%
  rename(
    IS01 = IS1,
    IS02 = IS2,
    IS03 = IS3,
    IS04 = IS4,
    IS05 = IS5,
    IS06 = IS6,
    IS07 = IS7,
    IS08 = IS8
  )

# clean and delete columns from Atlantic

eDNA_Actnow <- eDNA_Actnow %>%
  select(-c(PO01, PO02, PO03, PO04, PO05, PO06, PO07, PO08, PO09, PO10,
    PO11, PO12, PO13, PO14, PO15, WT03
  ))

## zero-out all reads ≤ 10 (including blanks)

no_modif <- c("id", "best_identity", "best_match", "order_name", "family_name", 
              "genus_name", "species_name", "scientific_name", "rank", 
              "sequence", "taxid")
change_cols <- setdiff(colnames(eDNA_Actnow), no_modif)

eDNA_Actnow[, change_cols] <- lapply(
  eDNA_Actnow[, change_cols, drop = FALSE],
  function(x) { if (is.numeric(x)) x[x <= 10] <- 0; x }
)

## blank subtraction
# subtract read counts from the associated sample columns
# per MPA_Identi for Filter_ID by PCR_ID

# subsets
metadata_Actnow <- metadata %>%
  filter(Source == "ACTNOW")
sample_meta <- metadata_Actnow %>%
  filter(str_detect(Filter_ID, "^Blank_"))
blank_meta <- metadata_Actnow %>%
  filter(str_detect(Filter_ID, "^Blank_"))

sample_pcrs <- sample_meta$PCR_ID
blank_pcrs  <- unique(blank_meta$PCR_ID)

# function to process one MPA
process_mpa <- function(mpa_id) {
  samples <- sample_meta %>%
    filter(MPA_Identi == mpa_id) %>%
    pull(PCR_ID)
  blanks <- blank_meta %>%
    filter(MPA_Identi == mpa_id) %>%
    pull(PCR_ID) %>%
    unique()
  if (length(samples) == 0 || length(blanks) == 0) return(NULL)
  blank_signal <- rowSums(eDNA_Actnow[, blanks, drop = FALSE], na.rm = TRUE) 
  # ^ extraction and pcr blanks summed together into one contamination vector 
  # (usually just one with reads)
  sample_mat <- eDNA_Actnow[, samples, drop = FALSE]
  sample_mat <- sweep(sample_mat, 1, blank_signal, FUN = "-")
  sample_mat[sample_mat < 0] <- 0
  eDNA_Actnow[, samples] <<- sample_mat
}

# run for all MPAs
mpas <- unique(metadata_Actnow$MPA_Identi)
walk(mpas, process_mpa)

# remove blank columns
eDNA_Actnow <- eDNA_Actnow %>%
  select(-all_of(blank_pcrs))

### write files

# dir.create("3_1_eDNA_initial_cleaned", showWarnings = FALSE)

write.csv(eDNA_Actnow, "3_1_eDNA_initial_cleaned/3_1_eDNA_Actnow.csv", row.names = FALSE)
write.csv(eDNA_Spygen, "3_1_eDNA_initial_cleaned/3_1_eDNA_Spygen.csv", row.names = FALSE)

### organise and match columns of both files

first_cols <- c(
  "id",
  "best_identity",
  "best_match",
  "taxid",
  "sequence",
  "class_name",
  "order_name",
  "family_name",
  "genus_name",
  "species_name",
  "scientific_name",
  "rank"
)

remaining_cols <- setdiff(names(eDNA_Actnow), first_cols)
eDNA_Actnow <- eDNA_Actnow[, c(first_cols, remaining_cols)]

cols_to_remove <- c(
  "count",
  "family",
  "genus",
  "order",
  "species",
  "id_status.db_gb264_and_custom2502_teleo",
  "rank_by_db.db_gb264_and_custom2502_teleo",
  "scientific_name_by_db.db_gb264_and_custom2502_teleo",
  "species_list.db_gb264_and_custom2502_teleo",
  "taxid_by_db.db_gb264_and_custom2502_teleo"
)
eDNA_Spygen <- eDNA_Spygen[, !(names(eDNA_Spygen) %in% cols_to_remove)]
remaining_cols <- setdiff(names(eDNA_Spygen), first_cols)
eDNA_Spygen <- eDNA_Spygen[, c(first_cols, remaining_cols)]

### merge eDNA data

# add source
eDNA_Actnow$id <- paste0("ACTNOW_", eDNA_Actnow$id)
eDNA_Spygen$id <- paste0("SPYGEN_", eDNA_Spygen$id)

# merge
key_cols <- names(eDNA_Actnow)[2:12]
eDNA_merged <- full_join(
  eDNA_Actnow,
  eDNA_Spygen,
  by = key_cols
)

# clean
eDNA_merged <- eDNA_merged %>%
  select(-id.x, -id.y)
eDNA_merged <- eDNA_merged %>%
  mutate(across(12:ncol(eDNA_merged), ~replace_na(.x, 0)))

eDNA_merged$family_name <- gsub("Cepolidae \\(in: bony fishes\\)", "Cepolidae", eDNA_merged$family_name)


eDNA_merged <- eDNA_merged[!duplicated(eDNA_merged), ] 

write.csv(eDNA_merged, "3_1_eDNA_initial_cleaned/3_1_eDNA_merged.csv", row.names = FALSE)
