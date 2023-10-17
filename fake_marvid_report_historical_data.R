
library(dplyr)

data <- read.csv(file.path(Data,"species_densities.csv"))

data <- data %>% rename(taxon=clean_taxonomy, density = density_n100m2) %>% dplyr::select(SampID, taxon, density)



sample_info <- read.csv(file.path(Data,"sample_info.csv"))

sample_info <- sample_info %>% rename(video.line=VL, section=SampID2, mean.latitude=y_coord, mean.longitude=x_coord) %>%
  dplyr::select(video.line, section, mean.latitude, mean.longitude, SampID)
