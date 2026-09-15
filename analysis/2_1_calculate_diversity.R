# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("vegan")
# install.packages("betapart")

library(dplyr)
library(tidyr)
library(vegan)
library(betapart)

eDNA_merged_MPA_PA <- read.csv("0_1_presence_absence/eDNA_merged_MPA_PA.csv", header = TRUE, sep = ",", dec = ".", check.names = FALSE)
combined_database <- read.csv("../data/1_species_database/1_4_additions_database_combined/1_4_combined_database_complete.csv", header = TRUE, sep = ",", dec = ".")

### filter for invasive species only 

invasive_eDNA_merged_MPA_PA <- eDNA_merged_MPA_PA %>%
  left_join(
    combined_database %>%
      select(species_name_fishbase, invasive, origin),
    by = "species_name_fishbase"
  ) %>%
  filter(invasive == 1)

### cut MPAs with no invasive species

invasive_eDNA_merged_MPA_PA <- invasive_eDNA_merged_MPA_PA %>%
  select(where(~ !all(. == 0, na.rm = TRUE)))

# ### gamma diversity (total) = number of invasive species
# 
# gamma_diversity <- invasive_eDNA_merged_MPA_PA %>%
#   summarise(
#     gamma_diversity = n_distinct(species_name_fishbase)
#   )
# 
# ### alpha diversity (per location) = richness per MPA
# 
# alpha_diversity <- invasive_eDNA_merged_MPA_PA %>%
#   pivot_longer(
#     cols = matches("^\\d+"),
#     names_to = "MPA_ID",
#     values_to = "presence"
#   ) %>%
#   filter(presence > 0) %>%
#   group_by(MPA_ID) %>%
#   summarise(
#     alpha_diversity = n_distinct(species_name_fishbase),
#     .groups = "drop"
#   )

### beta diversity (between locations)

# matrix

mpa_matrix <- invasive_eDNA_merged_MPA_PA %>%
  select(matches("^\\d+")) %>%
  t()
mpa_matrix <- as.matrix(mpa_matrix)
# dim(mpa_matrix)

# Jaccard index

beta_indices <- betapart::beta.pair(
  mpa_matrix,
  index.family = "jaccard"
)

jac_diss <- beta_indices$beta.jac
jac_turn <- beta_indices$beta.jtu
jac_nest <- beta_indices$beta.jne

# function 

get_beta_val <- function(Data, sites, names){
  res <- as.matrix(Data[[1]])
  rownames(res) <- colnames(res) <- sites
  res[upper.tri(res, diag = TRUE)] <- NA
  res_tot <- na.omit(
    cbind(
      expand.grid(dimnames(res)),
      Mean = as.vector(res)
    )
  )
  colnames(res_tot)[3] <- paste(names, "Mean", sep="_")
  res_tot
}

# output 

sites_list <- rownames(mpa_matrix)

Beta_jac <- get_beta_val(
  Data = list(jac_diss),
  sites = sites_list,
  names = "jac_diss"
)

Beta_turn <- get_beta_val(
  Data = list(jac_turn),
  sites = sites_list,
  names = "jac_turn"
)

Beta_nes <- get_beta_val(
  Data = list(jac_nest),
  sites = sites_list,
  names = "jac_nest"
)

average_beta_fd_indices <- Beta_jac %>% # different lengths cause some MPAs no invasive species
  left_join(
    Beta_turn %>% select(Var1, Var2, jac_turn_Mean),
    by = c("Var1", "Var2")
  ) %>%
  left_join(
    Beta_nes %>% select(Var1, Var2, jac_nest_Mean),
    by = c("Var1", "Var2")
  )

rownames(average_beta_fd_indices) <- paste(
  average_beta_fd_indices[,1],
  average_beta_fd_indices[,2],
  sep="_"
)

### overall mean and SD of beta diversity

overall_beta_stats <- data.frame(
  index = c("jac_diss", "jac_turn", "jac_nest"),
  mean = c(
    mean(Beta_jac$jac_diss_Mean, na.rm = TRUE),
    mean(Beta_turn$jac_turn_Mean, na.rm = TRUE),
    mean(Beta_nes$jac_nest_Mean, na.rm = TRUE)
  ),
  sd = c(
    sd(Beta_jac$jac_diss_Mean, na.rm = TRUE),
    sd(Beta_turn$jac_turn_Mean, na.rm = TRUE),
    sd(Beta_nes$jac_nest_Mean, na.rm = TRUE)
  )
)

overall_beta_stats

###

# dir.create("2_diversity", showWarnings = FALSE)

write.csv(average_beta_fd_indices , "2_diversity/2_1_beta_diversity.csv", row.names = FALSE)
write.csv(overall_beta_stats , "2_diversity/2_1_overall_beta_stats.csv", row.names = FALSE)

