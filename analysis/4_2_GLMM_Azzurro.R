# install.packages("tidyr")
# install.packages("dplyr")
# install.packages("glmmTMB")
# install.packages("broom.mixed")
# install.packages("purrr")
# install.packages("ggplot2")
# install.packages("RANN")
# install.packages("ncdf4")
# install.packages("gdistance")
# install.packages("sp")

library(tidyr)
library(dplyr)
library(glmmTMB)
library(broom.mixed)
library(purrr)
library(ggplot2)
library(RANN)
library(ncdf4)
library(raster)
library(gdistance)
library(sp)

azzurro_simple <- read.csv("3_comparison/3_azzurro_simple.csv", header = TRUE, sep = ",", dec = ".", check.names = FALSE)
threat_med <- readRDS("Threat_med_Coll.rds")
nc_sst <- nc_open("sst_mon_ltm_1991_2020.nc")
nc_sss <- nc_open("cmems_mod_glo_phy_anfc_0.083deg_PT1H-m_1787130902893.nc")
nc_grav <- nc_open("gravity4326.nc")
nc_circ <- nc_open("cmems_obs-mob_glo_phy-cur_my_0.25deg_P1M-m_1787830666253.nc")

###

azzurro_simple$eventDate <- as.Date(azzurro_simple$eventDate)
azzurro_simple$year <- format(azzurro_simple$eventDate, "%Y")

# keep only 2014 onward (checked with histogram)
azzurro_simple <- azzurro_simple %>%
  filter(
    !is.na(year),
    as.numeric(as.character(year)) >= 2014
  )

# grouping variables = rounding coordinates
azzurro_simple <- azzurro_simple %>%
  mutate(
    lat_group = round(decimalLatitude, 1),
    lon_group = round(decimalLongitude, 1),
    spatial_group = interaction(lat_group, lon_group)
  )

# invasive richness per spatial group
azzurro_simple <- azzurro_simple %>%
  group_by(spatial_group) %>%
  mutate(invasive_richness = n_distinct(species_name_fishbase)) %>%
  ungroup()

# cleaning of duplicates later before GLMM loop

### extrapolation 

## extract data coords

coords <- azzurro_simple %>%
  mutate(
    row_id = row_number()
  ) %>%
  dplyr::select(
    row_id,
    lat = decimalLatitude,
    lon = decimalLongitude
  )

## function

extrapolate_func <- function(
    lon,
    lat,
    lon_grid,
    lat_grid,
    mat,
    radius = 5
) {
  
  # nearest grid cell
  lon_i <- which.min(abs(lon_grid - lon))
  lat_i <- which.min(abs(lat_grid - lat))
  
  # local window
  lon_inds <- max(1, lon_i - radius):
    min(length(lon_grid), lon_i + radius)
  
  lat_inds <- max(1, lat_i - radius):
    min(length(lat_grid), lat_i + radius)
  
  local_mat <- mat[lon_inds, lat_inds]
  
  valid <- which(!is.na(local_mat), arr.ind = TRUE)
  
  if (nrow(valid) < 10) {
    return(NA_real_)
  }
  
  x <- lon_grid[lon_inds][valid[, 1]]
  y <- lat_grid[lat_inds][valid[, 2]]
  z <- local_mat[valid]
  
  # local quadratic surface
  fit <- tryCatch(
    lm(
      z ~ x + y + I(x^2) + I(y^2) + I(x*y)
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(NA_real_)
  }
  
  newdata <- data.frame(
    x = lon,
    y = lat
  )
  
  as.numeric(predict(fit, newdata = newdata))
}

### SSS data

# coordinates
lon_grid <- ncvar_get(nc_sss, "longitude")
lat_grid <- ncvar_get(nc_sss, "latitude")

# time
time_raw   <- ncvar_get(nc_sss, "time")
time_units <- ncatt_get(nc_sss, "time", "units")$value

origin <- sub("seconds since ", "", time_units)

file_date <- as.POSIXct(
  time_raw,
  origin = origin,
  tz = "UTC"
)

sss_full <- ncvar_get(nc_sss, "so")

# close NetCDF
nc_close(nc_sss)

# matrices
sss_mat <- sss_full

## function

# extrapolate

sss_data <- coords %>%
  rowwise() %>%
  mutate(
    sss = extrapolate_func(
      lon = lon,
      lat = lat,
      lon_grid = lon_grid,
      lat_grid = lat_grid,
      mat = sss_mat,
      radius = 5
    )
  ) %>%
  ungroup()

# add to azzurro data

azzurro_simple <- azzurro_simple %>%
  mutate(
    row_id = row_number()
  ) %>%
  left_join(
    sss_data %>%
      dplyr::select(row_id, sss),
    by = "row_id"
  ) %>%
  dplyr::select(-row_id)

### SST data 

lon_grid <- ncvar_get(nc_sst, "lon")   
lat_grid <- ncvar_get(nc_sst, "lat")

time_raw   <- ncvar_get(nc_sst, "time")
time_units <- ncatt_get(nc_sst, "time", "units")$value 
origin     <- sub(".*since ", "", time_units)
file_dates <- as.Date(time_raw, origin = origin)        

sst_full <- ncvar_get(nc_sst, "sst")

# close 
nc_close(nc_sst)

# matrices
sst_mat <- sst_full[, , 1]

# extrapolate

sst_data <- coords %>%
  rowwise() %>%
  mutate(
    sst = extrapolate_func(
      lon = lon,
      lat = lat,
      lon_grid = lon_grid,
      lat_grid = lat_grid,
      mat = sst_mat,
      radius = 5
    )
  ) %>%
  ungroup()

# add to azzurro data

azzurro_simple <- azzurro_simple %>%
  mutate(
    row_id = row_number()
  ) %>%
  left_join(
    sst_data %>%
      dplyr::select(row_id, sst),
    by = "row_id"
  ) %>%
  dplyr::select(-row_id)

### gravity data 

lon_grid <- ncvar_get(nc_grav, "lon")
lat_grid <- ncvar_get(nc_grav, "lat")
grav_full <- ncvar_get(nc_grav, "Band1")  
nc_close(nc_grav)

# matrices
grav_mat <- grav_full

# extrapolate

grav_data <- coords %>%
  rowwise() %>%
  mutate(
    grav = extrapolate_func(
      lon = lon,
      lat = lat,
      lon_grid = lon_grid,
      lat_grid = lat_grid,
      mat = grav_mat,
      radius = 5
    )
  ) %>%
  ungroup()

# add to azzurro data

azzurro_simple <- azzurro_simple %>%
  mutate(
    row_id = row_number()
  ) %>%
  left_join(
    grav_data %>%
      dplyr::select(row_id, grav),
    by = "row_id"
  ) %>%
  dplyr::select(-row_id)

### threat data

for (df_name in names(threat_med)) {
  
  threat_df <- threat_med[[df_name]]
  
  names(threat_df)[names(threat_df) == "X"] <- "decimalLongitude"
  names(threat_df)[names(threat_df) == "Y"] <- "decimalLatitude"
  
  threat_df <- threat_df %>%
    dplyr::select(decimalLongitude, decimalLatitude, var1.pred)
  
  threat_coords <- as.matrix(
    threat_df[, c("decimalLongitude", "decimalLatitude")]
  )
  
  sample_coords <- as.matrix(
    azzurro_simple[, c("decimalLongitude", "decimalLatitude")]
  )
  
  nn <- nn2(
    threat_coords,
    sample_coords,
    k = 1
  )
  
  azzurro_simple[[df_name]] <- 
    threat_df$var1.pred[nn$nn.idx[, 1]]
}

## clean for analysis

variables_to_remove <- c(
  "eventDate",
  "Oil spills",                                   
  "Risk of hypoxia",  
  "Invasive species",
  "Sea surface temperature increase",             
  "Coastal population density",  
  "Benthic structures oil rigs",  
  "Urban runoff",   
  "Organic pollution pesticides",                 
  "UV radiation", 
  "Nutrient input fertilizers",
  "Ocean acidification",
  "Commercial Shipping"
)

azzurro_simple <- azzurro_simple[, !(names(azzurro_simple) %in% variables_to_remove), drop = FALSE]

## mean of all fishing varibales

azzurro_simple$Fishing_Mean <- rowMeans(azzurro_simple[, c(
  "Fishing demersal non destructive high bycatch",
  "Fishing demersal non destructive low bycatch",
  "Fishing demersal destructive",
  "Fishing pelagic low bycatch",
  "Fishing pelagic high bycatch",
  "Artisanal fishing"
)], na.rm = TRUE)

azzurro_simple <- azzurro_simple[, !names(azzurro_simple) %in% c(
  "Fishing demersal non destructive high bycatch",
  "Fishing demersal non destructive low bycatch",
  "Fishing demersal destructive",
  "Fishing pelagic low bycatch",
  "Fishing pelagic high bycatch",
  "Artisanal fishing"
)]

### circulation data

cur_lon_grid <- ncvar_get(nc_circ, "longitude")
cur_lat_grid <- ncvar_get(nc_circ, "latitude")
u_stack <- drop(ncvar_get(nc_circ, "uo")) 
v_stack <- drop(ncvar_get(nc_circ, "vo"))
nc_close(nc_circ)

# uses scalar speed only (isotropic cost)
speed <- sqrt(u_stack^2 + v_stack^2)

# # sin / cos direction
# theta   <- atan2(v_stack, u_stack)
# dir_cos <- cos(theta)
# dir_sin <- sin(theta)

# cost distance from suez canal
suez_coord <- data.frame(lon = 32.32, lat = 31.27)

r_speed <- raster(t(speed[, ncol(speed):1]),
                  xmn = min(cur_lon_grid), xmx = max(cur_lon_grid),
                  ymn = min(cur_lat_grid), ymx = max(cur_lat_grid),
                  crs = "+proj=longlat +datum=WGS84")

conductance <- r_speed
vals <- values(conductance)
vals[vals <= 0 | is.na(vals)] <- 0.001
values(conductance) <- vals

tr <- transition(conductance, transitionFunction = mean, directions = 8)
tr <- geoCorrection(tr, type = "c")

suez_pt  <- SpatialPoints(suez_coord, proj4string = CRS("+proj=longlat +datum=WGS84"))
site_pts <- SpatialPoints(azzurro_simple[, c("decimalLongitude", "decimalLatitude")],
                          proj4string = CRS("+proj=longlat +datum=WGS84"))

azzurro_simple$circulation_cost_dist_suez <- as.numeric(costDistance(tr, suez_pt, site_pts))

### pre GLMM work

## clean up unrealistic extrapolations

azzurro_simple <- azzurro_simple %>%
  mutate(
    sst = na_if(sst, 0),
    sss = na_if(sss, 0),
    grav = na_if(grav, 0)
  )

## scale predictors

azzurro_final <- azzurro_simple %>%
  mutate(
    across(
      9:ncol(.),
      ~ if(is.numeric(.)) as.numeric(scale(.)) else .
    )
  )

predictors <- names(azzurro_final)[9:ncol(azzurro_final)]
predictors <- predictors[sapply(azzurro_final [predictors], is.numeric)]

## one row per spatial group 
azzurro_final <- azzurro_final %>%
  group_by(spatial_group) %>%
  summarise(
    invasive_richness = n_distinct(species_name_fishbase),
    lat_group = mean(lat_group, na.rm = TRUE),
    lon_group = mean(lon_group, na.rm = TRUE),
    across(all_of(predictors), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

### GLMM loop

uni_glmm_results_azzurro <- map_dfr(predictors, function(var){
  
  form <- as.formula(
    paste0("invasive_richness ~ `", var, "` + (1 | spatial_group)")
  )
  
  model <- glmmTMB(
    form,
    family = poisson,
    data = azzurro_final
  )
  
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(component == "cond", term != "(Intercept)") %>%
    mutate(
      predictor = var,
      AIC = AIC(model)
    ) %>%
    dplyr::select(predictor, estimate, std.error, statistic, p.value, AIC)
})

# rank
uni_glmm_results_azzurro <- uni_glmm_results_azzurro %>%
  arrange(p.value)

### coefficent plot

# significance labels
uni_glmm_results_azzurro <- uni_glmm_results_azzurro %>%
  mutate(
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

# add confidence intervals
uni_glmm_results_azzurro <- uni_glmm_results_azzurro %>%
  mutate(
    conf.low  = estimate - 1.96 * std.error,
    conf.high = estimate + 1.96 * std.error
  )

# rename predictors
uni_glmm_results_azzurro<- uni_glmm_results_azzurro %>%
  mutate(predictor = case_when(
    predictor == "sst" ~ "Mean Sea Surface Temperature (SST)",
    predictor == "sss" ~ "Mean Sea Surface Salinity (SSS)", 
    predictor == "Fishing_Mean" ~ "Fishing Gravity",
    predictor == "grav" ~ "Human Gravity",
    predictor == "circulation_cost_dist_suez" ~ "Circulation Cost-Distance to Red Sea",
    TRUE ~ predictor
  ))

# plot

plot_uni_glmm_azzurro <- ggplot(
  uni_glmm_results_azzurro,
  aes(
    y = reorder(predictor, estimate),
    x = estimate
  )
) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             colour = "grey50") +
  geom_errorbar(
    aes(
      xmin = conf.low,
      xmax = conf.high
    ),
    orientation = "y",
    width = 0.2,
    linewidth = 0.8
  ) +
  geom_point(size = 3) +
  geom_text(
    aes(
      x = max(conf.high, na.rm = TRUE) + 0.05,
      label = sig
    ),
    hjust = 0,
    size = 4,
    colour = "black"
  ) +
  labs(
    x ="Coefficient Estimate (Log Scale)",
    y = NULL
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        legend.position = "right")

print(plot_uni_glmm_azzurro)

ggsave(
  filename = "4_GLMM/4_2_azzurro_coefficent_plot_GLMM.pdf",
  plot = plot_uni_glmm_azzurro,
  device = "pdf",
  width = 10,
  height = 3,
  dpi = 1200
)

write.csv(uni_glmm_results_azzurro , "4_GLMM/4_2_azzurro_glmm_results_table.csv", row.names = FALSE)
