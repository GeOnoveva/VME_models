
library(dplyr)

data <- read.csv(file.path(Data,"species_densities.csv"))

data <- data %>% rename(taxon=clean_taxonomy, density = density_n100m2) %>% select(SampID, taxon, density)



sample_info <- read.csv(file.path(Data,"sample_info.csv"))

sample_info <- sample_info %>% rename(video.line=VL, section=SampID2, mean.latitude=y_coord, mean.longitude=x_coord) %>%
  select(video.line, section, mean.latitude, mean.longitude, SampID)



## fake spatial pixels data frame (rather, incomplete)

i<-1

predpix <- raster(paste0(ftp_filepath,
                         paste(unlist(strsplit(files[i], split="all", fixed=TRUE))[2], sep=""))
)
names(predpix) <- substr(paste(unlist(strsplit(files[i], split="all", fixed=TRUE))[2], sep=""),
                         2,
                         nchar(paste(unlist(strsplit(files[i], split="all", fixed=TRUE))[2], sep=""))-4)
predpix <- aggregate(predpix,4)

predpix <- as(predpix, "SpatialPixelsDataFrame")


i<-2

tmpraster <- raster(paste0(ftp_filepath,
                           paste(unlist(strsplit(files[i], split="all", fixed=TRUE))[2], sep=""))
)

tmpraster <- aggregate(tmpraster,4)


predpix[[i]]<- values(tmpraster)[predpix@grid.index]

names(predpix@data)[[2]]<-substr(paste(unlist(strsplit(files[i], split="all", fixed=TRUE))[2], sep=""),
                                 2,
                                 nchar(paste(unlist(strsplit(files[i], split="all", fixed=TRUE))[2], sep=""))-4)
