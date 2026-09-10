# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("ggplot2")

library(dplyr)
library(tidyr)
library(ggplot2)

invasive_richness_MPA_red_sea <- read.csv("1_richness/1_1_invasive_richness_MPA_red_sea.csv", header = TRUE, sep = ",", dec = ".")
invasive_richness_MPA_atlantic <- read.csv("1_richness/1_1_invasive_richness_MPA_atlantic.csv", header = TRUE, sep = ",", dec = ".")
sample_info <- read.csv("../data/2_Mediterranean_map/sample_locations_info.csv", header = TRUE, sep = ",", dec = ".")

### summary of needed values

summary <- sample_info %>%
  dplyr::select(
    MPA_Identi,
    decimalLatitude,
    decimalLongitude
  )

plot_df <- summary %>%
  left_join(
    invasive_richness_MPA_atlantic,
    by = c("MPA_Identi" = "MPA_ID")
  ) %>%
  left_join(
    invasive_richness_MPA_red_sea,
    by = c("MPA_Identi" = "MPA_ID")
  ) %>%
  mutate(
    invasive_richness_MPA_atlantic = replace_na(invasive_richness_MPA_atlantic, 0),
    invasive_richness_MPA_red_sea = replace_na(invasive_richness_MPA_red_sea, 0)
  ) %>%
  pivot_longer(
    cols = c(
      invasive_richness_MPA_atlantic,
      invasive_richness_MPA_red_sea
    ),
    names_to = "origin",
    values_to = "richness"
  )

### plot

ggplot(plot_df, aes(
  x = decimalLongitude,
  y = richness,
  colour = origin
)) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    span = 0.5,       # higher = smoother
    linewidth = 1.5
  ) + 
  geom_point(
    size = 2.5,
    alpha = 0.5
  ) +
  scale_colour_manual(
    values = c(
      "invasive_richness_MPA_atlantic" = "#4687a1",
      "invasive_richness_MPA_red_sea" = "#d8613e"
    ),
    labels = c(
      "invasive_richness_MPA_atlantic" = "Atlantic Ocean",
      "invasive_richness_MPA_red_sea" = "Red Sea"
    )
  ) +
  scale_x_continuous(
    labels = function(x) {
      ifelse(x < 0,
             paste0(abs(x), "°W"),
             ifelse(x == 0,
                    "0°",
                    paste0(x, "°E")))
    }
  ) +
  labs(
    x = "Longitude",
    y = "Non-Native Species Richness"
  ) +
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)
        # , panel.grid = element_blank()
  )

# dir.create("1_richness", recursive = TRUE)

ggsave(
  filename = "1_richness/1_3_richness_scatterplot.pdf",  
  plot = last_plot(),                         
  device = "pdf",                             
  width = 4,                                  
  height = 4,                                 
  dpi = 1200                                   
)