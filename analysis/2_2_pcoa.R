# install.packages("ade4")
# install.packages("sf")
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("rnaturalearth")
# install.packages("scales")

library(ade4)
library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(scales)

beta_diversity <- read.csv("2_diversity/2_1_beta_diversity.csv", header = TRUE, sep = ",", dec = ".")
combined_database <- read.csv("../data/1_species_database/1_2_combined_database_complete.csv", header = TRUE, sep = ",", dec = ".")
metadata <- read.csv("../data/3_eDNA_cleaning/METADATA_MERGED_edited.csv", header = TRUE, sep = ",", dec = ".")
sample_locations <- read.csv("../data/2_Mediterranean_map/sample_locations_info.csv")

# move one malta and one marseille location slightly so that it is clear on the map
sample_locations <- sample_locations %>%
  mutate(
    decimalLongitude = if_else(MPA_Identi == 40001, decimalLongitude + 0.1, decimalLongitude),
    decimalLatitude = if_else(MPA_Identi == 40001, decimalLatitude - 0.1, decimalLatitude),
    decimalLongitude = if_else(MPA_Identi == 90009, decimalLongitude + 0.1, decimalLongitude),
    decimalLatitude = if_else(MPA_Identi == 90009, decimalLatitude - 0.1, decimalLatitude)
  )

# sample locations for MPAs with no invasive species
sample_locations_empty <- sample_locations %>%
  filter(
    MPA_Identi %in% c(50001, 50002, 50003, 80001, 80002, 90005, 90007, 90009, 70006, 70002, 70003)
  )

### pivot to square distance matrices 

all_sites <- sort(unique(c(as.character(beta_diversity$Var1),
                           as.character(beta_diversity$Var2))))
n <- length(all_sites)

make_dist_matrix <- function(long_df, value_col) {
  mat <- matrix(0, nrow = n, ncol = n,
                dimnames = list(all_sites, all_sites))
  for (i in seq_len(nrow(long_df))) {
    r <- as.character(long_df$Var1[i])
    c <- as.character(long_df$Var2[i])
    v <- long_df[[value_col]][i]
    if (!is.na(v)) {
      mat[r, c] <- v
      mat[c, r] <- v
    } else {
      mat[r, c] <- 1
      mat[c, r] <- 1
    }
  }
  as.dist(mat)
}

dist_jac  <- make_dist_matrix(beta_diversity, "jac_diss_Mean")  # total Jaccard
dist_turn <- make_dist_matrix(beta_diversity, "jac_turn_Mean")  # turnover
dist_nest <- make_dist_matrix(beta_diversity, "jac_nest_Mean")  # nestedness

### PCoA

PCoA_jac <- dudi.pco(quasieuclid(dist_jac), scannf = FALSE, nf = 8)
Eigen_jac <- PCoA_jac$eig / sum(PCoA_jac$eig) * 100

PCoA_turn <- dudi.pco(quasieuclid(dist_turn), scannf = FALSE, nf = 8)
Eigen_turn <- PCoA_turn$eig / sum(PCoA_turn$eig) * 100

PCoA_nest <- dudi.pco(quasieuclid(dist_nest), scannf = FALSE, nf = 8)
Eigen_nest <- PCoA_nest$eig / sum(PCoA_nest$eig) * 100

# extract coordinates
dat_jac <- as.data.frame(PCoA_jac$li[, 1:2])
dat_turn <- as.data.frame(PCoA_turn$li[, 1:2])
dat_nest <- as.data.frame(PCoA_nest$li[, 1:2])

# change colnames
colnames(dat_jac) <- c("PCoA1", "PCoA2")
colnames(dat_turn) <- c("PCoA1", "PCoA2")
colnames(dat_nest) <- c("PCoA1", "PCoA2")

# add site ids
dat_jac$MPA_Identi <- rownames(dat_jac)
dat_turn$MPA_Identi <- rownames(dat_turn)
dat_nest$MPA_Identi <- rownames(dat_nest)

# add region info
metadata$MPA_Identi <- as.character(metadata$MPA_Identi)
dat_jac$MPA_Identi <- as.character(dat_jac$MPA_Identi)
dat_turn$MPA_Identi <- as.character(dat_turn$MPA_Identi)
dat_nest$MPA_Identi <- as.character(dat_nest$MPA_Identi)

metadata_region <- metadata %>%
  select(MPA_Identi, Country) %>%
  distinct()

dat_jac <- left_join(dat_jac,
                     metadata_region,
                     by = "MPA_Identi")

dat_turn <- left_join(dat_turn,
                      metadata_region,
                      by = "MPA_Identi")

dat_nest <- left_join(dat_nest,
                      metadata_region,
                      by = "MPA_Identi")

### plot pcoa

dat_jac$colour_value_jac <- atan2(dat_jac$PCoA2, dat_jac$PCoA1)
dat_turn$colour_value_turn <- atan2(dat_turn$PCoA2, dat_turn$PCoA1)
dat_nest$colour_value_nest <- atan2(dat_nest$PCoA2, dat_nest$PCoA1)

rainbow_pal <- hue_pal()(12) 

## add ecoregion info

ecoregions <- sample_locations %>%
  mutate(MPA_Identi = as.character(MPA_Identi)) %>%
  select(MPA_Identi, Ecoregion) %>%
  distinct()

dat_jac  <- left_join(dat_jac,  ecoregions, by = "MPA_Identi")
dat_turn <- left_join(dat_turn, ecoregions, by = "MPA_Identi")
dat_nest <- left_join(dat_nest, ecoregions, by = "MPA_Identi")

## jac

plot_jac <- ggplot(dat_jac, aes(
    x = PCoA1, 
    y = PCoA2, 
    colour = colour_value_jac,
    fill = colour_value_jac,
    shape = Ecoregion
  )) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_point(size = 6) +
  scale_shape_manual(
    values = c(
      "Levantine Sea" = 24,
      "Aegean Sea" = 23,
      "Ionian Sea" = 22,
      "Western Mediterranean" = 25
    )
  ) +
  scale_colour_gradientn(
    colours = rainbow_pal
  ) +
  scale_fill_gradientn(
    colours = rainbow_pal
  ) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.margin = margin(5, 15, 5, 5)
  ) +
  labs(
    x = paste0("PCoA1 (", round(Eigen_jac[1],1), "%)"),
    y = paste0("PCoA2 (", round(Eigen_jac[2],1), "%)")
  )
print(plot_jac)

## turn

plot_turn <- ggplot(dat_turn, aes(x = PCoA1, y = PCoA2, 
                                colour = colour_value_turn,
                                fill = colour_value_turn,
                                shape = Ecoregion
)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_point(size = 6) +  
  scale_shape_manual(
    values = c(
      "Levantine Sea" = 24,
      "Aegean Sea" = 23,
      "Ionian Sea" = 22,
      "Western Mediterranean" = 25
    )
  ) +
  scale_colour_gradientn(
    colours = rainbow_pal
  ) +
  scale_fill_gradientn(
    colours = rainbow_pal
  ) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.margin = margin(5, 15, 5, 5)
  ) +
  labs(
    x = paste0("PCoA1 (", round(Eigen_turn[1],1), "%)"),
    y = paste0("PCoA2 (", round(Eigen_turn[2],1), "%)")
  )
print(plot_turn)

### plot maps

world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_make_valid()

sample_locations$MPA_Identi <- as.character(sample_locations$MPA_Identi)

## jac

dat_jac_map <- left_join(
  dat_jac,
  sample_locations %>%
    select(MPA_Identi, decimalLongitude, decimalLatitude) %>%
    distinct(),
  by = "MPA_Identi"
)

map_jac <- ggplot() +
  
  geom_sf(
    data = world,
    fill = "#F5F0EB",
    color = "#A08A6F",
    linewidth = 0.2
  ) +
  
  geom_point(
    data = sample_locations_empty,
    aes(
      x = decimalLongitude,
      y = decimalLatitude
    ),
    shape = 16,     
    fill = NA,
    colour = "black",
    stroke = 0.6,
    size = 3.5,
    alpha = 0.4
  ) +
  
  geom_point(
    data = dat_jac_map,
    aes(
      x = decimalLongitude,
      y = decimalLatitude,
      colour = colour_value_jac,
      fill = colour_value_jac,
      shape = Ecoregion
    ),
    size = 4
  ) +
  
  scale_shape_manual(
    values = c(
      "Levantine Sea" = 24,
      "Aegean Sea" = 23,
      "Ionian Sea" = 22,
      "Western Mediterranean" = 25
    )
  ) +
  
  scale_colour_gradientn(
    colours = rainbow_pal,
    name = "PCoA Position"
  ) +
  scale_fill_gradientn(
    colours = rainbow_pal,
    guide = "none"
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
  theme(    legend.position = "right",
            legend.text = element_text(size = 11),
            legend.title = element_text(size = 12),
            legend.key = element_rect(fill = "white", colour = "white"),
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 13),
            panel.grid = element_blank(),
            panel.background = element_rect(fill = "#e7f5fb", colour = NA)
  )

print(map_jac)

## turn

dat_turn_map <- left_join(
  dat_turn,
  sample_locations %>%
    select(MPA_Identi, decimalLongitude, decimalLatitude) %>%
    distinct(),
  by = "MPA_Identi"
)

map_turn <- ggplot() +
  
  geom_sf(
    data = world,
    fill = "#F5F0EB",
    color = "#A08A6F",
    linewidth = 0.2
  ) +
  
  geom_point(
    data = sample_locations_empty,
    aes(
      x = decimalLongitude,
      y = decimalLatitude
    ),
    shape = 16,     
    fill = NA,
    colour = "black",
    stroke = 0.6,
    size = 3.5,
    alpha = 0.4
  ) +
  
  geom_point(
    data = dat_turn_map,
    aes(
      x = decimalLongitude,
      y = decimalLatitude,
      colour = colour_value_turn,
      fill = colour_value_turn,
      shape = Ecoregion
    ),
    size = 4
  ) +
  
  scale_shape_manual(
    values = c(
      "Levantine Sea" = 24,
      "Aegean Sea" = 23,
      "Ionian Sea" = 22,
      "Western Mediterranean" = 25
    )
  ) +
  
  scale_colour_gradientn(
    colours = rainbow_pal,
    name = "PCoA Position"
  ) +
  scale_fill_gradientn(
    colours = rainbow_pal,
    guide = "none"
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
  theme(    legend.position = "right",
            legend.text = element_text(size = 11),
            legend.title = element_text(size = 12),
            legend.key = element_rect(fill = "white", colour = "white"),
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 13),
            panel.grid = element_blank(),
            panel.background = element_rect(fill = "#e7f5fb", colour = NA)
  )

map_turn <- map_turn +
  guides(
    colour = guide_colourbar(order = 1),
    shape = guide_legend(order = 2)
  )

print(map_turn)

### save plots

ggsave(
  "2_diversity//2_2_PCoA_Jaccard_plot_jac.pdf",
  plot = plot_jac,
  width = 6,   
  height = 5,  
)

ggsave(
  "2_diversity//2_2_PCoA_Jaccard_map_jac.pdf",
  plot = map_jac,
  width = 9,   
  height = 5,  
)

ggsave(
  "2_diversity//2_2_PCoA_Jaccard_plot_turn.pdf",
  plot = plot_turn,
  width = 6,   
  height = 5,  
)

ggsave(
  "2_diversity//2_2_PCoA_Jaccard_map_turn.pdf",
  plot = map_turn,
  width = 9,   
  height = 5,  
)
