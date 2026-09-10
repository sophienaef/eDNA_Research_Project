# install.packages("tidyr")
# install.packages("dplyr")
# install.packages("sf")
# install.packages("geosphere")
# install.packages("ggplot2")
# install.packages("rnaturalearth")
# install.packages("scales")
# install.packages("patchwork")
# install.packages("png")
# install.packages("grid")

library(tidyr)
library(dplyr)
library(sf)
library(geosphere)
library(ggplot2)
library(rnaturalearth)
library(scales)
library(patchwork)
library(png)
library(grid)

eDNA_merged_Filter_ID_PA <- read.csv("0_1_presence_absence/eDNA_merged_Filter_ID_PA.csv", header = TRUE, sep = ",", dec = ".", check.names = FALSE)
metadata <- read.csv("../data/3_eDNA_cleaning/METADATA_MERGED_edited.csv", header = TRUE, sep = ",", dec = ".")
combined_database <- read.csv("../data/1_species_database/1_2_combined_database_complete.csv", header = TRUE, sep = ",", dec = ".")
ORMEF_occurrence_filtered <- read.csv("../data/1_species_database//1_1_ORMEF_occurrence_filtered.csv", header = TRUE, sep = ",", dec = ".")

Alepes_djedaba_silhouette <- readPNG("silhouettes/Alepes_djedaba_silhouette.png")
Callionymus_filamentosus_silhouette <- readPNG("silhouettes/Callionymus_filamentosus_silhouette.png")
Upeneus_pori_silhouette <- readPNG("silhouettes/Upeneus_pori_silhouette.png")

### filter for invasive species only 

invasive_eDNA_merged_Filter_ID_PA <- eDNA_merged_Filter_ID_PA %>%
  left_join(
    combined_database %>%
      select(species_name_fishbase, invasive, origin),
    by = "species_name_fishbase"
  ) %>%
  filter(invasive == 1)

### re order df

# reshape into long format
presence_long <- invasive_eDNA_merged_Filter_ID_PA %>%
  pivot_longer(
    cols = starts_with(c("DAR", "SPY")),
    names_to = "Filter_ID",
    values_to = "presence"
  ) %>%
  filter(presence == 1)

metadata <- metadata %>%
  mutate(Filter_ID = gsub("-", ".", Filter_ID))

# join with metadata to add coordinates
result_df <- presence_long %>%
  left_join(
    metadata %>%
      select(Filter_ID, MPA_Identi, decimalLongitude, decimalLatitude),
    by = "Filter_ID"
  ) %>%
  select(
    species_name_fishbase,
    Filter_ID,
    MPA_Identi,
    decimalLongitude,
    decimalLatitude
  )

### add species_name_fishbase

ORMEF_occurrence_filtered <- merge(
  ORMEF_occurrence_filtered,
  combined_database[, c("species_name", "species_name_fishbase")],
  by.x = "scientificName",
  by.y = "species_name",
  all.x = TRUE  
)

## replace missing

ORMEF_occurrence_filtered$species_name_fishbase[is.na(ORMEF_occurrence_filtered$species_name_fishbase)] <-
  ORMEF_occurrence_filtered$scientificName[is.na(ORMEF_occurrence_filtered$species_name_fishbase)]

ORMEF_occurrence_filtered$species_name_fishbase[ORMEF_occurrence_filtered$species_name_fishbase == "Silhouetta aegyptia"] <-
  "Silhouettea aegyptia"

### simplify df

found_simple <- result_df [, c("species_name_fishbase", 
                               "decimalLatitude",
                               "decimalLongitude", 
                               "MPA_Identi")] 

azzurro_simple <- ORMEF_occurrence_filtered [, c("species_name_fishbase",
                                                 "decimalLatitude",
                                                 "decimalLongitude")] 

### save database

azzurro_simple_export <- ORMEF_occurrence_filtered [, c("species_name_fishbase",
                                                 "eventDate",
                                           "decimalLatitude",
                                           "decimalLongitude")] 

write.csv(azzurro_simple_export,
          "3_comparison/3_azzurro_simple.csv",
          row.names = FALSE)

### compare 

matches_list <- list()

# loop 200km

for(i in seq_len(nrow(found_simple))) {
  sp <- found_simple$species_name_fishbase[i]
  lat <- found_simple$decimalLatitude[i]
  lon <- found_simple$decimalLongitude[i]
  
  inv_sp <- azzurro_simple %>%
    filter(species_name_fishbase == sp)
  
  if(nrow(inv_sp) == 0) next
 
  d <- distHaversine(
    matrix(c(lon, lat), ncol = 2),
    as.matrix(inv_sp[, c("decimalLongitude", "decimalLatitude")])
  )
  
  nearby <- inv_sp[d <= 200000, ]
  
  if(nrow(nearby) > 0) {
    nearby$distance_km <- d[d <= 200000] / 1000
    nearby$found_MPA <- found_simple$MPA_Identi[i]
    nearby$found_latitude <- lat
    nearby$found_longitude <- lon
    
    matches_list[[i]] <- nearby
  }
}

matches_200km <- bind_rows(matches_list)
matches_200km <- matches_200km %>%
  mutate(
    species_name_fishbase = forcats::fct_reorder(
      species_name_fishbase,
      distance_km,
      .fun = median
    )
  )

### box plot for each species

boxplot <- ggplot(
  matches_200km,
  aes(
    x = species_name_fishbase,
    y = distance_km,
    fill = species_name_fishbase
  )
) +
  geom_jitter(
    aes(colour = species_name_fishbase),
    width = 0.15,
    size = 0.5,
    alpha = 1
  ) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.4
  ) +
  labs(
    x = "Species",
    y = "Distance (km)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.position = "none",
    plot.margin = margin(5.5, 5.5, 5.5, 30)
  )

print(boxplot)

# dir.create("3_comparison", showWarnings = FALSE)
ggsave("3_comparison/3_boxplot_distance_by_species.pdf", plot = boxplot, width = 12, height = 8, dpi = 600)

### maps

# set colours for species

species_order_median <- levels(matches_200km$species_name_fishbase)

species_cols <- setNames(
  scales::hue_pal()(length(species_order_median)),
  species_order_median
)

unique_species <- sort(unique(as.character(matches_200km$species_name_fishbase)))

plot_list <- list()

for (species in unique_species) {
  species_data <- matches_200km %>%
    filter(species_name_fishbase == species)
  
  found_data <- species_data %>%
    select(found_latitude, found_longitude) %>%
    distinct() %>%
    rename(Latitude = found_latitude, Longitude = found_longitude) %>%
    mutate(Type = "Found")
  
  azzurro_data <- species_data %>%
    select(decimalLatitude, decimalLongitude) %>%
    distinct() %>%
    rename(Latitude = decimalLatitude, Longitude = decimalLongitude) %>%
    mutate(Type = "Azzurro")
  
  combined_data <- bind_rows(found_data, azzurro_data)
  
  lat_min <- min(combined_data$Latitude) - 5
  lat_max <- max(combined_data$Latitude) + 5
  lon_min <- min(combined_data$Longitude) - 5
  lon_max <- max(combined_data$Longitude) + 5
  
  world <- ne_countries(scale = "medium", returnclass = "sf") %>%
    st_crop(xmin = lon_min, ymin = lat_min, xmax = lon_max, ymax = lat_max)

  map <- ggplot() +
    geom_sf(
      data = world,
      fill = "#F5F0EB",
      color = "#A08A6F",
      linewidth = 0.2
    ) +
    
    geom_point(
      data = filter(combined_data, Type == "Azzurro"),
      aes(x = Longitude, y = Latitude),
      color = species_cols[species],
      size = 5,
      alpha = 0.5
    ) +
    
    geom_point(
      data = filter(combined_data, Type == "Found"),
      aes(x = Longitude, y = Latitude),
      color = "black",
      size = 5
    
    ) +
    coord_sf(
      xlim = c(lon_min + 1.5, lon_max - 1.5),  
      ylim = c(lat_min + 1.5, lat_max - 1.5), 
      expand = FALSE
    ) +
    labs(
      title = paste(species),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_bw() +
    theme(axis.text = element_text(size = 17),
          axis.title = element_text(size = 19),
          plot.title = element_text(size = 19),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "#e7f5fb", colour = NA))
  
  plot_list[[species]] <- map
  
}

chunk_size <- 6
plot_chunks <- split(plot_list, ceiling(seq_along(plot_list)/chunk_size))

pdf(
  "3_comparison/supporting_information/combined_maps.pdf",
  width = 20,  
  height = 24, 
  onefile = FALSE
)

for (i in seq_along(plot_chunks)) {
  
  combined_plot <- wrap_plots(plot_chunks[[i]], ncol = 2) +
    plot_layout(heights = rep(1, 3))
  
  pdf(
    paste0(
      "3_comparison/supporting_information/3_combined_maps_page_",
      i,
      ".pdf"
    ),
    width = 20,
    height = 24
  )
  
  print(combined_plot)
  dev.off()
}

# filtered 

filtered_matches_200km <- matches_200km %>%
     filter(species_name_fishbase %in% c("Alepes djedaba", "Callionymus filamentosus", "Upeneus pori"))

filtered_unique_species <- sort(
  unique(as.character(filtered_matches_200km$species_name_fishbase))
)

species_settings <- list(
  "Alepes djedaba" = list(
    axis_text = 21,
    axis_title = 24,
    title = 24,
    point_found = 6,
    point_azzurro = 6,
    set_colour = "#DF8A00"
  ),
  
  "Callionymus filamentosus" = list(
    axis_text = 21,
    axis_title = 23,
    title = 23,
    point_found = 6,
    point_azzurro = 6,
    set_colour = "#A2A400"
  ),
  
  "Upeneus pori" = list(
    axis_text = 15,
    axis_title = 16,
    title = 16,
    point_found = 5,
    point_azzurro = 5,
    set_colour = "#00C0B8"
  )
)

for (species in filtered_unique_species) {
  
  species_data <- matches_200km %>%
    filter(species_name_fishbase == species)
  
  settings <- species_settings[[species]]
  
  found_data <- species_data %>%
    select(found_latitude, found_longitude) %>%
    distinct() %>%
    rename(Latitude = found_latitude, Longitude = found_longitude) %>%
    mutate(Type = "Found")
  
  azzurro_data <- species_data %>%
    select(decimalLatitude, decimalLongitude) %>%
    distinct() %>%
    rename(Latitude = decimalLatitude, Longitude = decimalLongitude) %>%
    mutate(Type = "Azzurro")
  
  combined_data <- bind_rows(found_data, azzurro_data)
  
  lat_min <- min(combined_data$Latitude) - 5
  lat_max <- max(combined_data$Latitude) + 5
  lon_min <- min(combined_data$Longitude) - 5
  lon_max <- max(combined_data$Longitude) + 5
  
  world <- ne_countries(scale = "medium", returnclass = "sf") %>%
    st_crop(xmin = lon_min, ymin = lat_min, xmax = lon_max, ymax = lat_max)
  
  silhouette <- switch(
    species,
    "Alepes djedaba" = Alepes_djedaba_silhouette,
    "Callionymus filamentosus" = Callionymus_filamentosus_silhouette,
    "Upeneus pori" = Upeneus_pori_silhouette
  )
  
  silhouette_grob <- rasterGrob(
    silhouette,
    interpolate = TRUE
  )
  
  map <- ggplot() +
    geom_sf(
      data = world,
      fill = "#F5F0EB",
      color = "#A08A6F",
      linewidth = 0.2
    ) +
    
    annotation_custom(
      silhouette_grob,
      xmin = lon_min + 1.5,
      xmax = lon_max - 1.5,
      ymin = lat_min + 1.5,
      ymax = lat_max - 1.5
    ) +
    
    geom_point(
      data = filter(combined_data, Type == "Azzurro"),
      aes(x = Longitude, y = Latitude),
      color = settings$set_colour,
      size = settings$point_azzurro,
      alpha = 0.5
    ) +
    
    geom_point(
      data = filter(combined_data, Type == "Found"),
      aes(x = Longitude, y = Latitude),
      color = "black",
      size = settings$point_found
    ) +
    coord_sf(
      xlim = c(lon_min + 1.5, lon_max - 1.5),  
      ylim = c(lat_min + 1.5, lat_max - 1.5), 
      expand = FALSE
    ) +
    labs(
      title = paste(species),
      x = "Longitude",
      y = "Latitude"
    ) +
    scale_x_continuous(
      breaks = seq(
        ceiling(lon_min + 1.5),
        floor(lon_max - 1.5),
        by = ifelse(species == "Upeneus pori", 5, 2)
      )
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = settings$axis_text),
      axis.title = element_text(size = settings$axis_title),
      plot.title = element_text(size = settings$title),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "#e7f5fb", colour = NA)
    )
  
  ggsave(
    filename = paste0("3_comparison//3_map_", species, "_silhouette.pdf"),
    plot = map,
    width = 10,
    height = 8,
    dpi = 600
  )

}
