# analysis

## 0_1_presence_absence_conversion.R
converts raw eDNA data into presence-absence matrices

## 0_2_invasive_species_table.R
list of all non-native species found
produces Supporting Table 1

## 1_1_calculate_richness.R
calculations for further richness analysis

## 1_2_map_richness.R
produces Figure 1b

## 1_3_scatterplot_richness.R
produces Figure 1c

## 1_4_scatterplot_reads.R
produces Figure 1d

## 2_1_calculate_diversity.R
calculations for further diversity analysis
produces Supporting Table 2

## 2_2_pcoa.R
PCoA
produces Figure 2

## 3_azzurro_comparison.R
eDNA with database radius comparison
produces Figure 3 and Supporting Figures 6 to 10

## 4_1_GLMM_eDNA.R
GLMM of eDNA data
produces Figure 4b and Supporting Table 4

## 4_2_GLMM_Azzurro.R
GLMM of non-native species database (Azzurro et al. 2026)
produces Figure 4a and Supporting Table 3

## 5_read_plots.R
produces Supporting Figures 1 to 5

## inputs

### silhouettes directory
silhouettes for Figure 3
used in 3_azzurro_comparison.R

### cmems_mod_glo_phy_anfc_0.083deg_PT1H-m_1787130902893.nc (zipped)
SSS data
used in 4_2_GLMM_Azzurro.R

### cmems_obs-mob_glo_phy-cur_my_0.25deg_P1M-m_1787830666253.nc
circulation data
used in 4_1_GLMM_eDNA.R and 4_2_GLMM_Azzurro.R

### gravity4326.nc (zipped)
human gravity data 
used in 4_1_GLMM_eDNA.R and 4_2_GLMM_Azzurro.R

### sst_mon_ltm_1991_2020.nc (zipped)
SST data
used in 4_2_GLMM_Azzurro.R

### Threat_med_Coll.rds
threats / anthropogenic varibles
used in 4_1_GLMM_eDNA.R and 4_2_GLMM_Azzurro.R

## outputs
0_1_presence_absence directory
0_2_invasive_species_table directory
1_richness directory
2_diversity directory
3_comparison directory
4_GLMM directory
5_read_plots directory
