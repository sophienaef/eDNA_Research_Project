# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("ggplot2")

library(dplyr)
library(tidyr)
library(ggplot2)

eDNA_merged_MPA <- read.csv("../data/3_eDNA_cleaning/3_3_eDNA_region_cleaned/3_3_eDNA_merged_MPA.csv", header = TRUE, sep = ",", dec = ".", check.names = FALSE)
combined_database <- read.csv("../data/1_species_database/1_2_combined_database_complete.csv", header = TRUE, sep = ",", dec = ".")
sample_info <- read.csv("../data/2_Mediterranean_map/sample_locations_info.csv", header = TRUE, sep = ",", dec = ".")

### filter for invasive species only 

invasive_eDNA_merged_MPA <- eDNA_merged_MPA %>%
  left_join(
    combined_database %>%
      dplyr::select(species_name_fishbase, invasive, origin),
    by = "species_name_fishbase"
  ) %>%
  filter(invasive == 1)

### add up reads per MPA regardless of species only origin

read_counts <- invasive_eDNA_merged_MPA %>%
  dplyr::group_by(origin) %>%
  dplyr::summarise(
    dplyr::across(
      `10001`:`90009`,
      ~sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

read_counts_long <- read_counts %>%
  pivot_longer(
    cols = -origin,
    names_to = "MPA_Identi",
    values_to = "reads"
  ) %>%
  mutate(MPA_Identi = as.integer(MPA_Identi))

### add coordinates

plot_df_reads <- read_counts_long %>%
  left_join(
    sample_info %>% 
      dplyr::select(MPA_Identi, decimalLongitude) %>%
      distinct(),
    by = "MPA_Identi"
  )

### log reads

plot_df_reads <- plot_df_reads %>%
  mutate(
    reads_sqrt = sqrt(reads)
  )

### plot

ggplot(plot_df_reads, aes(
  x = decimalLongitude,
  y = reads_sqrt,
  colour = origin
)) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    span = 0.5,
    linewidth = 1.5
  ) +
  geom_point(
    size = 2.5,
    alpha = 0.5
  ) +
  scale_colour_manual(
    values = c(
      "A" = "#4687a1",
      "RS" = "#d8613e"
    ),
    labels = c(
      "A" = "Atlantic Ocean",
      "RS" = "Red Sea"
    )
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
  labs(
    x = "Longitude",
    y = "Square Root of Number of Reads \nof Non-Native Species",
    colour = "Origin of \nNon-Native \nSpecies"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

ggsave(
  filename = "1_richness/1_4_reads_scatterplot.pdf",
  plot = last_plot(),
  device = "pdf",
  width = 6,
  height = 4,
  dpi = 1200
)