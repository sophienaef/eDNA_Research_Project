# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("ggplot2")
# install.packages("patchwork")

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

eDNA_merged_MPA <- read.csv("../data/3_eDNA_cleaning/3_3_eDNA_region_cleaned/3_3_eDNA_merged_MPA.csv", header = TRUE, sep = ",", dec = ".", check.names = FALSE)
combined_database <- read.csv("../data/1_species_database/1_2_combined_database_complete.csv", header = TRUE, sep = ",", dec = ".")
metadata <- read.csv("../data/3_eDNA_cleaning/METADATA_MERGED_edited.csv", header = TRUE, sep = ",", dec = ".") 

### filter for invasive species only 

invasive_eDNA_merged_MPA <- eDNA_merged_MPA  %>%
  left_join(
    combined_database %>%
      select(species_name_fishbase, invasive, origin),
    by = "species_name_fishbase"
  ) %>%
  filter(invasive == 1)

### simplify data

columns_to_remove <- c("scientific_name", 
                       "rank", 
                       "species_name", 
                       "genus_name", 
                       "genus_name_fishbase",
                       "family_name", 
                       "invasive", 
                       "origin")

invasive_eDNA_merged_MPA <- invasive_eDNA_merged_MPA %>%
  select(-all_of(columns_to_remove))

### convert data to long fomat

invasive_long <- invasive_eDNA_merged_MPA %>%
  pivot_longer(
    cols = -species_name_fishbase,
    names_to = "MPA_Identi",
    values_to = "reads"
  ) %>%
  mutate(MPA_Identi = as.character(MPA_Identi))

# add coordinates

metadata_clean <- metadata %>%
  mutate(MPA_Identi = as.character(MPA_Identi)) %>%
  distinct(MPA_Identi, .keep_all = TRUE)

invasive_long <- invasive_long %>%
  left_join(
    metadata_clean %>%
      select(MPA_Identi, decimalLongitude, decimalLatitude),
    by = "MPA_Identi"
  )

### plots

# dir.create("5_read_plots", showWarnings = FALSE)

species_list <- sort(unique(invasive_long$species_name_fishbase))

## bar plots - one per species

plot_list <- list()

for (sp in species_list) {
  
  plot_data <- invasive_long %>%
    filter(species_name_fishbase == sp) %>%
    arrange(desc(reads))

  p <- ggplot(plot_data, aes(
    x = decimalLongitude,
    y = reads
  )) +
    geom_col(
      width = 0.15
    ) +
    labs(
      title = sp,
      x = "Longitude",
      y = "Number of eDNA Reads"
    ) +
    scale_x_continuous(
      labels = function(x) {
        ifelse(
          x < 0,
          paste0(abs(x), "°W"),
          ifelse(
            x == 0,
            "0°",
            paste0(x, "°E")
          )
        )
      }
    ) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 13),
      axis.title = element_text(size = 15),
      plot.title = element_text(size = 15)
    )
  
  plot_list[[sp]] <- p
}

### combine plots into pages

chunk_size <- 6

plot_chunks <- split(
  plot_list,
  ceiling(seq_along(plot_list) / chunk_size)
)


for (i in seq_along(plot_chunks)) {
  
  combined_plot <- wrap_plots(
    plot_chunks[[i]],
    ncol = 2
  ) +
    plot_layout(
      heights = rep(1, 3)
    )
  
  ggsave(
    filename = paste0(
      "5_read_plots/5_combined_read_plots_page_",
      i,
      ".pdf"
    ),
    plot = combined_plot,
    width = 14,
    height = 15,
    dpi = 600
  )
}
