# install.packages("sf")
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("rnaturalearth")
# install.packages("png")
# install.packages("grid")

library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(png)
library(grid)

sample_locations <- read.csv("sample_locations_info.csv")
arrow_layer <- readPNG("arrow_layer_clear.png")

# load world map
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_make_valid()

meow <- st_read("meow/meow_ecos.shp") %>%
  st_make_valid()

bathy <- st_read("bathy/med_etopo2_proj.shp") %>%
  st_make_valid() %>%
  st_transform(st_crs(world))

#

# move one malta and one marseille location slightly so that it is clear on the map
sample_locations <- sample_locations %>%
  mutate(
    decimalLongitude = if_else(MPA_Identi == 40001, decimalLongitude + 0.1, decimalLongitude),
    decimalLatitude = if_else(MPA_Identi == 40001, decimalLatitude - 0.1, decimalLatitude),
    decimalLongitude = if_else(MPA_Identi == 90009, decimalLongitude + 0.1, decimalLongitude),
    decimalLatitude = if_else(MPA_Identi == 90009, decimalLatitude - 0.1, decimalLatitude)
  )

med_names <- c(
  "Alboran Sea",
  "Western Mediterranean",
  "Adriatic Sea",
  "Ionian Sea",
  "Aegean Sea",
  "Levantine Sea",
  "Tunisian Plateau/Gulf of Sidra"
)

med <- meow %>%
  filter(ECOREGION %in% med_names) %>%
  st_transform(st_crs(world))

# sampling locations

unique_transects <- sample_locations %>%
  distinct(
    decimalLatitude,
    decimalLongitude,
    MPA_Identi
  ) %>%
  filter(!is.na(decimalLatitude))

transects_sf <- st_as_sf(
  unique_transects,
  coords = c("decimalLongitude", "decimalLatitude"),
  crs = 4326
) %>%
  st_transform(st_crs(world))

# load ecoregions
meow <- st_read("meow/meow_ecos.shp") %>%
  st_make_valid()

med_names <- c(
  "Alboran Sea",
  "Western Mediterranean",
  "Adriatic Sea",
  "Ionian Sea",
  "Aegean Sea",
  "Levantine Sea",
  "Tunisian Plateau/Gulf of Sidra"
)

med <- meow %>%
  filter(ECOREGION %in% med_names) %>%
  st_transform(st_crs(world))

# load depths / bathymetry

bathy_sea <- bathy %>%
  filter(DN < 0)
bathy_coast <- bathy %>%
  filter(DN < 0 & DN >= -200)
bathy_coast_outline <- bathy_coast %>%
  st_union() %>%
  st_as_sf()

# default label positions
med_labels_ecoregions <- st_point_on_surface(med)

# manually adjusted labels
custom_labels <- data.frame(
  ECOREGION = c(
    "Alboran \nSea",
    "Western \nMediterranean",
    "Adriatic Sea",
    "Ionian Sea",
    "Aegean Sea",
    "Gulf of Sidra"
  ),
  # x = c(-1.5, 6.1, 17.6, 19, 25.8, 17.5),
  # y = c(36.3, 39.7, 42, 37, 36, 33.5)
  x = c(-0.95, 6.2, 15.25, 19, 25.8, 17.5),
  y = c(36.55, 40.7, 43.1, 37.8, 36, 33.4)
) %>%
  st_as_sf(coords = c("x", "y"), crs = st_crs(world))


original_names <- c(
  "Alboran Sea",
  "Western Mediterranean",
  "Adriatic Sea",
  "Ionian Sea",
  "Aegean Sea",
  "Tunisian Plateau/Gulf of Sidra"
)

combined_ecoregion_labels <- bind_rows(
  med_labels_ecoregions %>% 
    filter(!ECOREGION %in% original_names),
  custom_labels
)

# sea and strait labels

label_data <- data.frame(
  name = c(
    "Red Sea",
    "Atlantic \nOcean",
    "Strait of \nGibraltar",
    "Suez Canal",
    "Strait of \nSicily"
  ),
  
  x = c(
    34.6, -7.0, -5.35, 31.2, 12.8
  ),
  
  y = c(
    27.3, 35.4, 36.65, 30.4, 36.85
  )
)

# dahsed lines for straits

connection_lines <- data.frame(
  name = c("Strait of Gibraltar", "Suez Canal", "Strait of Sicily"),
  x_start = c(-6.25, 32.54,12.2),
  y_start = c(35.8, 29.75, 36.75),
  x_end   = c(-5.0, 32.5, 11.2),
  y_end   = c(36.05, 31.35, 37.8)
)

# plot

ggplot() +
  
  # depths / bathymetry
  geom_sf(
    data = bathy_coast,
    fill = "#f6fbfe",
    colour = NA
  ) +
  
  # ecoregions
  geom_sf(
    data = med,
    fill = NA,
    color = "#006B50", #568A5A
    linewidth = 0.18
  ) +
  
  # world and countries
  geom_sf(
    data = world,
    fill = "#F5F0EB",
    color = "#A08A6F",
    linewidth = 0.2
  ) +   
  
  # samples
  geom_sf(
    data = transects_sf,
    color = "black", #73b5c7 #FF8C00 
    size = 1.6,
    shape = 19
  ) +
  
  geom_segment(
    data = connection_lines,
    aes(
      x = x_start,
      y = y_start,
      xend = x_end,
      yend = y_end
    ),
    linetype = "dashed",
    linewidth = 0.25,
    color = "#1a5f7a" #006AA0
  ) +
  
  # sea and strait labels
  geom_text(
    data = label_data,
    aes(
      x = x,
      y = y,
      label = name
    ),
    size = 1.8,
    color = "#1a5f7a",
    vjust = 0
  ) +
  
  # ecoregion labels
  geom_sf_text(
    data = combined_ecoregion_labels,
    aes(
      label = ECOREGION,
      geometry = geometry
    ),
    size = 1.8,
    color = "#006B50"
  ) +
  
  annotation_custom(
    rasterGrob(
      arrow_layer,
      interpolate = TRUE
    ),
    xmin = -7.35,
    xmax = 35.7,
    ymin = 26.1,
    ymax = 46.75
  ) +
  
  coord_sf(
    xlim = c(-8, 37),
    ylim = c(27, 47),
    expand = FALSE
  ) +
  
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "#e7f5fb", colour = NA)
  )

# dir.create("2_output_maps", recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = "2_output_maps/2_1_labels_combined_map.pdf",
  plot = last_plot(),
  device = "pdf",
  width = 8,
  height = 5,
  dpi = 1200
)