# install.packages("tidyr")
# install.packages("dplyr")
# install.packages("glmmTMB")
# install.packages("broom.mixed")
# install.packages("purrr")
# install.packages("ggplot2")
# install.packages("corrplot")
# install.packages("car")
# install.packages("factoextra")
# install.packages("RANN")
# install.packages("ncdf4")
# install.packages("tidyverse")
# install.packages("sf")
# install.packages("raster")
# install.packages("gdistance")
# install.packages("sp")

library(tidyr)
library(dplyr)
library(glmmTMB)
library(broom.mixed)
library(purrr)
library(ggplot2)
library(corrplot)
library(car)
library(factoextra)
library(RANN)
library(ncdf4)
library(tidyverse)
library(sf)
library(raster)
library(gdistance)
library(sp)

invasive_richness_Transect_ID <- read.csv("1_richness/1_1_invasive_richness_Transect_ID.csv", header = TRUE, sep = ",", dec = ".", check.names = FALSE)
metadata <- read.csv("../data/3_eDNA_cleaning/METADATA_MERGED_edited.csv", header = TRUE, sep = ",", dec = ".")
threat_med <- readRDS("Threat_med_Coll.rds")
nc_sst <- nc_open("sst_mon_ltm_1991_2020.nc")
nc_sss <- nc_open("cmems_mod_glo_phy_anfc_0.083deg_PT1H-m_1787130902893.nc")
nc_grav <- nc_open("gravity4326.nc")
nc_circ <- nc_open("cmems_obs-mob_glo_phy-cur_my_0.25deg_P1M-m_1787830666253.nc")

### filter metadata

short_metadata <- metadata[, c(
  "Transect_ID", "MPA_Identi",
  "decimalLatitude", "decimalLongitude"
)]

# remove blank rows
short_metadata <- short_metadata %>%
  filter(!if_any(everything(), is.na))

# remove duplicates info per Transect_ID
short_metadata <- short_metadata %>%
  distinct()
# n_distinct(short_metadata$Transect_ID)

short_metadata <- short_metadata %>%
  filter(Transect_ID != "90007.2.2")

analysis_data <- short_metadata

### extrapolation 

## extract data coords

coords <- analysis_data %>%
  dplyr::select(
    Transect_ID,
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

# add to data

analysis_data <- analysis_data %>%
  left_join(
    sss_data %>%
      dplyr::select(Transect_ID, sss),
    by = "Transect_ID"
  )

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

# add to data

analysis_data <- analysis_data %>%
  left_join(
    sst_data %>%
      dplyr::select(Transect_ID, sst),
    by = "Transect_ID"
  )

### gravity data (extrapolation)

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

## add to data

analysis_data <- analysis_data %>%
  left_join(
    grav_data %>%
      dplyr::select(Transect_ID, grav),
    by = "Transect_ID"
  )

### add threat info

short_metadata$decimalLatitude <- round(short_metadata$decimalLatitude, 2)
short_metadata$decimalLongitude <- round(short_metadata$decimalLongitude, 2)

# # split threat_med into dfs
# list2env(threat_med, envir = .GlobalEnv)

enriched_metadata <- short_metadata

for (df_name in names(threat_med)) {
  
  threat_df <- threat_med[[df_name]]
  
  names(threat_df)[names(threat_df) == "X"] <- "decimalLongitude"
  names(threat_df)[names(threat_df) == "Y"] <- "decimalLatitude"
  
  threat_df <- threat_df %>%
    dplyr::select(decimalLongitude, decimalLatitude, var1.pred)
  
  threat_coords <- as.matrix(threat_df[, c("decimalLongitude", "decimalLatitude")])
  sample_coords <- as.matrix(enriched_metadata[, c("decimalLongitude", "decimalLatitude")])
  
  nn <- nn2(threat_coords, sample_coords, k = 1)
  
  enriched_metadata[[df_name]] <- threat_df$var1.pred[nn$nn.idx]
}

### merge data for analysis

analysis_data <- analysis_data %>%
  left_join(
    invasive_richness_Transect_ID %>%
      rename(invasive_richness = invasive_richness_Transect_ID),
    by = "Transect_ID"
  ) %>%
  mutate(
    invasive_richness = replace_na(invasive_richness, 0)
  ) %>%
  relocate(invasive_richness, .after = MPA_Identi)

# remove duplicate transects
analysis_data <- analysis_data %>%
  distinct(Transect_ID, .keep_all = TRUE)

# add enriched

enriched_metadata <- enriched_metadata %>%
  dplyr::select(-MPA_Identi, -decimalLatitude, -decimalLongitude)

analysis_data <- analysis_data %>%
  left_join(
    enriched_metadata,
    by = "Transect_ID"
  )

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
site_pts <- SpatialPoints(analysis_data[, c("decimalLongitude", "decimalLatitude")],
                          proj4string = CRS("+proj=longlat +datum=WGS84"))

analysis_data$circulation_cost_dist_suez <- as.numeric(costDistance(tr, suez_pt, site_pts))

#####

### variable reduction

## remove additional variables 

variables_to_remove <- c(
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
  "Commercial Shipping",
  "Gravity"
)

analysis_data <- analysis_data[, !(names(analysis_data) %in% variables_to_remove), drop = FALSE]

## mean of all fishing variables

analysis_data$Fishing_Mean <- rowMeans(analysis_data[, c(
  "Fishing demersal non destructive high bycatch",
  "Fishing demersal non destructive low bycatch",
  "Fishing demersal destructive",
  "Fishing pelagic low bycatch",
  "Fishing pelagic high bycatch",
  "Artisanal fishing"
)], na.rm = TRUE)

analysis_data <- analysis_data[, !names(analysis_data) %in% c(
  "Fishing demersal non destructive high bycatch",
  "Fishing demersal non destructive low bycatch",
  "Fishing demersal destructive",
  "Fishing pelagic low bycatch",
  "Fishing pelagic high bycatch",
  "Artisanal fishing"
)]

analysis_scaled <- analysis_data %>%
  mutate(
    across(
      6:ncol(.),
      ~ if(is.numeric(.)) as.numeric(scale(.)) else .
    )
  )

predictors <- names(analysis_scaled)[6:ncol(analysis_scaled)]
predictors <- predictors[sapply(analysis_scaled[predictors], is.numeric)]

### GLMM loop

uni_glmm_results <- map_dfr(predictors, function(var){
  
  form <- as.formula(
    paste0("invasive_richness ~ `", var, "` + (1 | MPA_Identi)")
  )
  
  model <- glmmTMB(
    form,
    family = poisson,
    data =  analysis_scaled
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
uni_glmm_results <- uni_glmm_results %>%
  arrange(p.value)

# significance labels
uni_glmm_results <- uni_glmm_results %>%
  mutate(
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

# add confidence intervals
uni_glmm_results <- uni_glmm_results %>%
  mutate(
    conf.low  = estimate - 1.96 * std.error,
    conf.high = estimate + 1.96 * std.error
  )

# rename predictors
uni_glmm_results <- uni_glmm_results %>%
  mutate(predictor = case_when(
    predictor == "sst" ~ "Mean Sea Surface Temperature (SST)",
    predictor == "sss" ~ "Mean Sea Surface Salinity (SSS)", 
    predictor == "Fishing_Mean" ~ "Fishing Gravity",
    predictor == "grav" ~ "Human Gravity",
    predictor == "circulation_cost_dist_suez" ~ "Circulation Cost-Distance to Red Sea",
    TRUE ~ predictor
  ))

uni_glmm_results <- uni_glmm_results %>%
  mutate(predictor = factor(predictor, levels = c("Circulation Cost-Distance to Red Sea",
                                                  "Fishing Gravity",
                                                  "Human Gravity",
                                                  "Mean Sea Surface Salinity (SSS)",
                                                  "Mean Sea Surface Temperature (SST)"
  ))) %>%
  arrange(predictor)

# plot

plot_uni_glmm <- ggplot(
  uni_glmm_results,
  aes(
    y = predictor,
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
      x = max(conf.high, na.rm = TRUE) + 0.3,
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

print(plot_uni_glmm)

ggsave(
  filename = "4_GLMM/4_1_eDNA_coefficent_plot_GLMM.pdf",
  plot = plot_uni_glmm,
  device = "pdf",
  width = 10,
  height = 3,
  dpi = 1200
)

# write csv

uni_glmm_results <- uni_glmm_results %>%
  mutate(across(where(is.numeric), ~ round(.x, 5))) %>%
  arrange(factor(predictor, levels = c(
    "Mean Sea Surface Temperature (SST)",
    "Mean Sea Surface Salinity (SSS)",
    "Human Gravity",
    "Fishing Gravity",
    "Circulation Cost-Distance to Red Sea"
  )))

write.csv(uni_glmm_results , "4_GLMM/4_1_eDNA_glmm_results_table.csv", row.names = FALSE,
          quote = FALSE)
