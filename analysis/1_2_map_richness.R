# install.packages("sf")
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("rnaturalearth")
# install.packages("tidyr")
# install.packages("gstat")
# install.packages("raster")

library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(tidyr)
library(gstat)
library(raster)

invasive_richness_MPA <- read.csv("1_richness/1_1_invasive_richness_MPA.csv", header = TRUE, sep = ",", dec = ".")
sample_info <- read.csv("../data/2_Mediterranean_map/sample_locations_info.csv", header = TRUE, sep = ",", dec = ".")

### summary of needed values

summary <- sample_info %>%
  dplyr::select(
    MPA_Identi,
    decimalLatitude,
    decimalLongitude
  )
summary <- summary %>%
  left_join(
    invasive_richness_MPA,
    by = c("MPA_Identi" = "MPA_ID")
  )
summary <- summary %>%
  mutate(
    invasive_richness_MPA = replace_na(invasive_richness_MPA, 0)
  )

summary_sf <- st_as_sf(
  summary,
  coords = c("decimalLongitude", "decimalLatitude"),
  crs = 4326
)

# convert richness points to Spatial
summary_sp <- as(summary_sf, "Spatial")

# create prediction grid
r <- raster(
  extent(-8, 37, 27, 47),
  res = 0.1,
  crs = CRS("+proj=longlat +datum=WGS84")
)

grid <- rasterToPoints(r, spatial = TRUE)

# create world and cut for med

world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_make_valid()

# interpolate richness

idw_richness <- idw(
  invasive_richness_MPA ~ 1,
  summary_sp,
  newdata = grid,
  idp = 2 # how conservative interpolation is
)

idw_df <- as.data.frame(idw_richness)
names(idw_df)[1:3] <- c("x", "y", "richness")

idw_sf <- st_as_sf(
  idw_df,
  coords = c("x", "y"),
  crs = 4326
)

# approximate mediterranean boundary
med_polygon <- st_polygon(list(
  matrix(
    c(
      -6, 36,    # Gibraltar
      0, 43.5,   # Spain/France coast
      10, 44.5,  # Ligurian/North Med
      19, 48,    # Adriatic
      27, 40,    # Turkey
      36, 37,    # Levant
      35, 30,    # Egypt
      25, 31,    # south Levant
      15, 25,    # Libya
      5, 28,     # Tunisia
      -5, 35,    # Algeria/Morocco
      -6, 36     # close
    ),
    ncol = 2,
    byrow = TRUE
  )
)) %>%
  st_sfc(crs = 4326) %>%
  st_as_sf()

idw_med <- st_filter(
  idw_sf,
  med_polygon,
  .predicate = st_within
)

idw_med_df <- cbind(
  st_coordinates(idw_med),
  st_drop_geometry(idw_med)
)

names(idw_med_df)[1:2] <- c("x", "y")

# plot

ggplot() +
  
  # interpolated richness
  geom_raster(
    data = idw_med_df,
    aes(
      x = x,
      y = y,
      fill = richness
    )
  ) +
  
  geom_sf(
    data = world,
    fill = "#F5F0EB",
    color = "#A08A6F",
    linewidth = 0.2
  ) +
  
  scale_fill_gradientn(
    name = "Non-Native \nSpecies \nRichness",
    colours = c(
      "#1a5f7a",
      "#8cc3d0",
      "#d9e8e8",
      "#f5d4a8",
      "#e88a4e",
      "#c8372d"
    ),
    values = scales::rescale(c(
      0, 1, 2, 5, 12, 19
    )),
    limits = c(0, 19),
    breaks = c(0, 2, 5, 10, 19),
    labels = c("0", "2", "5", "10", "19"),
    guide = guide_colorbar(
      barheight = unit(5, "cm"),
      barwidth = unit(0.6, "cm"),
      ticks = TRUE
    )
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
  theme(legend.position = "right",
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "#e7f5fb", colour = NA)
  )

# dir.create("1_richness", recursive = TRUE)

ggsave(
  filename = "1_richness/1_2_richness_map.pdf",  
  plot = last_plot(),                         
  device = "pdf",                             
  width = 8,                                  
  height = 5,                                 
  dpi = 1200                                   
)
