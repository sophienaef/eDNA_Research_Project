# install.packages("sf")
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("rnaturalearth")
# install.packages("ncdf4")

library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(ncdf4)

nc_sss_data <- nc_open("cmems_mod_glo_phy_anfc_0.083deg_PT1H-m_1787130902893.nc")

salinity_data <- ncvar_get(nc_sss_data, "so")  
lon <- ncvar_get(nc_sss_data, "longitude")
lat <- ncvar_get(nc_sss_data, "latitude")

nc_close(nc_sss_data)

# crop data
lon_min <- -8
lon_max <- 37
lat_min <- 27
lat_max <- 47
lon_idx <- which(lon >= lon_min & lon <= lon_max)
lat_idx <- which(lat >= lat_min & lat <= lat_max)

salinity_med <- salinity_data[
  lon_idx,
  lat_idx
]
lon_med <- lon[lon_idx]
lat_med <- lat[lat_idx]

sss_df <- expand.grid(
  lon = lon_med,
  lat = lat_med
)
sss_df$sss <- as.vector(salinity_med)
sss_df <- na.omit(sss_df)

world <- ne_countries(scale = "medium", returnclass = "sf")

# plot

ggplot() +
  
  geom_tile(data = sss_df, aes(x = lon, y = lat, fill = sss)) +
  
  geom_sf(
    data = world,
    fill = "#F5F0EB",
    color = "#A08A6F",
    linewidth = 0.2
  ) +  
  scale_fill_gradientn(
    name = "SSS (PSU)",
    colors = c(
      "#104257",
      "#1a5f7a",
      "#73b5c7",
      "#a5d0d9",
      "#d9e8e8",
      "#f5d4a8",
      "#e88a4e",
      "#c8372d"
    ),
    
    values = c(
      0.00, # 0
      0.76, # 35
      0.84, # 37.5
      0.92, # 40
      1.00  # 42.5
    ),
    
    limits = c(0, 42.5),
    
    breaks = c(
      0,
      10,
      20,
      30,
      34,
      36,
      38,
      40,
      42
    ),
    
    labels = c(
      "0",
      "10",
      "20",
      "30",
      "34",
      "36",
      "38",
      "40",
      "42"
    ),
    
    na.value = "#104257",
    
    guide = guide_colorbar(
      barheight = unit(7.5, "cm"),
      barwidth = unit(0.6, "cm"),
      ticks = TRUE
    )
  )+
  coord_sf(
    xlim = c(lon_min, lon_max),
    ylim = c(lat_min, lat_max),
    expand = FALSE
  ) +
  
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "SSS (PSU)"
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
  filename = "2_output_maps/2_3_SSS_map.pdf",
  plot = last_plot(),
  device = "pdf",
  width = 8,
  height = 5,
  dpi = 600
)
