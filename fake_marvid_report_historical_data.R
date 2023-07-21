
library(dplyr)

old_data <- read.csv("U:\\Mareano\\VIDEOLAB\\VIDEO DATA\\200m_scale_species_by_sample\\Data_Delivery_2024\\species_densities.csv")

old_data <- old_data %>% rename(taxon=clean_taxonomy, density = density_n100m2) %>% select(SampID, taxon, density)

data <- old_data


old_sample_info <- read.csv("U:\\Mareano\\VIDEOLAB\\VIDEO DATA\\200m_scale_species_by_sample\\Data_Delivery_2023\\inputs\\sample_info.csv")

old_sample_info <- old_sample_info %>% rename(video.line=VL, section=SampID2, mean.latitude=y_coord, mean.longitude=x_coord) %>%
  select(video.line, section, mean.latitude, mean.longitude, SampID)

sample_info <- old_sample_info

# IMPORTANT NOTE: I may have to re-generate sample_info and species_densities (from the make matrix generalized notebook!!)

rm(old_data, old_sample_info)
