# install.packages("sf")
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("rnaturalearth")

library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)

sst_2024 <- read.csv("avg_sst_med_2024.csv")

world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_make_valid()

# plot

ggplot() +
  
  geom_tile(data = sst_2024, aes(x = lon, y = lat, fill = avg_sst)) +
  
  geom_sf(
    data = world,
    fill = "#F5F0EB",
    color = "#A08A6F",
    linewidth = 0.2
  ) +  
  scale_fill_gradientn(
      name = "SST (°C)",
      colors = c(
        "#1a5f7a",
        "#73b5c7",
        "#a5d0d9",
        "#d9e8e8",
        "#f5d4a8",
        "#e88a4e",
        "#c8372d"
      ),
    #   colors = c(
    #     "#0000FF",
    #     "#00FFFF",
    #     "#00FF00",
    #     "#FFFF00",
    #     "#FF7F00",
    #     "#FF0000"
    #   ),
    # values = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    limits = c(13.34199, 27.05494), 
  ) +
  
  coord_sf(
    xlim = c(-8, 37),
    ylim = c(27, 47),
    expand = FALSE
  ) +
  
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "SST (°C)"
  ) +
  
  theme_bw() +
  theme(legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "#e7f5fb", colour = NA)
        )

# dir.create("2_output_maps", recursive = TRUE)

ggsave(
  filename = "2_output_maps/2_2_SST_map.pdf",  
  plot = last_plot(),                         
  device = "pdf",                             
  width = 8,                                  
  height = 5,                                 
  dpi = 1200                                   
)

