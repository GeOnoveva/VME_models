

## Remember to wipe the environment when running from the top!

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc()

# Paths:

## Path to data folder on network drive:
Data <- "G:/R_projects/VH_MarVid/Data"
#Data <- "~/R_projects/VH_MarVid/Data" # when running from Server

## Path to HI OneDrive folder, where Taxonary normally resides (only accessible when working from desktop RStudio)
#SharedData <- "C:\\Users\\genoveva\\Havforskningsinstituttet\\Video analysis lab - Reference List\\"


## Path to ftp site, where environmental layers (already cropped and aligned) reside
ftp_url <- "ftp://ftp2.ngu.no/toGeno/EnvLayersNiN2022/all"
ftp_filepath <- paste0(strsplit(ftp_url, "\\//")[[1]][1], "//havforsk:Havf8776@", 
                   strsplit(ftp_url, "\\//")[[1]][2])#insert the credentials into the url

## Path to folder where to save intermediate objects created along the way and needed in subsequent steps (you need to create this folder if it doesn't exist)
savepath <- "G:/R_projects/VH_MarVid/Intermediate_objects"
#savepath <- "~/R_projects/VH_MarVid/Intermediate_objects" # when running from Server


## Final products to be shared:
outpath <- "G:/R_projects/VH_MarVid/Results"
#outpath <- "~/R_projects/VH_MarVid/Results" # when running from Server


# Functions:

## function to generate sample ID

makesampid <- function(x){
  require(dplyr)
  x <- x %>% mutate(SampID=paste(video.line, section, sep = "_"))
}

## function to project to utm

utmize <- function(x){
  require(rgdal)
  proj4string(x)=CRS("+init=epsg:4326")
  result <- spTransform(x, CRS("+init=epsg:32633")) 
  return(result)
}

## functions to rescale between 5 and 260

rang <- function(x)                    
{
  rng <- rep(NA_real_, length(x))
  minx <- min(x, na.rm = T)
  maxx <- max(x, na.rm = T)
  for(i in seq(length(x))) {
    rng[i] <- (((x[i] - minx) / (maxx - minx)) * 255) + 5
  }
  return(rng)
}

rescl5_260 <- function(df)
{
  ## Copying metadata of raster stack
  rs5_260 <- df
  
  ## Looping over every raster in stack
  for(i in seq(dim(rs5_260)[2])) {
    var.i <- df[,i]
    scaled.i <- rang(var.i)
    rs5_260[,i] <- scaled.i
    #colnames(rs5_260)[i] <- colnames(var.i)
  }
  rs5_260
}

## function to remove from v observations where response in the low class, then split the rest randomly, all while keeping lines together

stratisplit <- function (x, n, keep, p, video.line) {
  x<-v
  n<-3
  keep=0.25
  p=0.3
  video.line=video.line # lookup video line by matching v to sample info on SampID
  #require(BAMMtools, groupdata2, tidyverse)
  
  # stratification by response
  densclas <- cut(v$vmeind_density, breaks = getJenksBreaks(v$vmeind_density, n+1),# 3 strata (jenks breaks)
                  include.lowest = TRUE,
                  labels = FALSE) 
  
  # first, find out which video lines have most sections in the low dens class: defined as having max 1 section in a higher class:
  d1 <- data.frame(cbind(densclas,video.line))
  d1$value <- 1
  d2 <- pivot_wider(d1, names_from = "densclas", values_fn = sum, values_fill = 0)
  colnames(d2)[2:4] <- c("low", "interm", "high")
  d2 <- d2 %>% mutate(d2, tot = rowSums(dplyr::select(.,low, interm, high)),
                      max = low/tot)
  d3 <- d2 %>% filter(max>0.5)
  
  # grab a small random subset of low dens lines (defined by the keep parameter)
  d3 <- partition(d3, keep, list_out = FALSE)
  
  # reserve partition = 2 (will count towards the final training set)
  v <- data.frame(cbind(video.line,v))
  v_core <- v %>% filter(!(video.line%in%d3$video.line[which(d3$.partitions==2)])) # video lines with medium and high densities of vme indicators
  desiredprop <- 1-p
  neededprop <- desiredprop-((dim(v)[1]-dim(v_core)[1])/100)
  
  # split the core normally
  v_split <- partition(v_core, p=1-neededprop, id_col = 'video.line', list_out = FALSE)
  
  # v_test = v_split test (partition = 2 is test)
  v_test = v_split %>% filter(.partitions==2)
  
  # v_train = vsplit train + those reserved from before (partition = 1 is train)
  d4 <- v %>% filter(video.line%in%d3$video.line[which(d3$.partitions==2)])
  v_train = v_split %>% filter(.partitions==1) %>%
    dplyr::select(-which(colnames(.)==".partitions")) %>% rbind(.,d4)
  return(list(v_test, v_train))
}

stratisplit2 <- function(x, p){
  x<-v
  p=0.2
  # stratification by response
  densclas <- cut(v$mean_dens, breaks = getJenksBreaks(v$vmeind_density, n+1),# 3 strata (jenks breaks)
                  include.lowest = TRUE,
                  labels = FALSE)
  v$video.line <- as.factor(sub("_.*", "", v$SampID))
  core <- v %>% filter(densclas%in%c(2,3))
  fringe <- v %>% filter(densclas==1)
  
  d1 <- partition(core, p=0.2, id_col= "video.line", list_out = FALSE)
  d2 <- partition(fringe, p=0.5, id_col= "video.line", list_out = FALSE)
  
  v_train <- rbind(filter(d1,.partitions==2),filter(d2,.partitions==2))
  v_test <- rbind(filter(d1,.partitions==1),filter(d2,.partitions==1))
  return(list(v_test, v_train))
}

## function to polygonize the high density predicted areas

threspol <- function(thresvalue, raster, min.area, max.hole){
  #require(SpatialEco, raster, smoothr, units, sf)
  rcl <- matrix(data = c(0,thresvalue,thresvalue,max(values(raster), na.rm =TRUE),0,1), nrow=2)
  rec <- reclassify(raster, rcl, include.lowest = TRUE)
  pol <- rasterToPolygons(rec, fun=function(x){x==1}, dissolve = TRUE)
  #pol <- as(spatialEco::explode(pol), "Spatial")
  pol<-sf::st_cast(as(pol, "sf"))
  pol<- as(pol, "Spatial")
  pol$area_sqkm <- area(pol) / 1000000
  pol <- pol[pol$area_sqkm>min.area,]
  area_thresh <- units::set_units(max.hole, km^2)
  pol <- smoothr::fill_holes(pol, threshold = area_thresh)
  pol <- smoothr::smooth(pol, method = "ksmooth")
  return(pol)
}

## function to derive buffer distances for a set of points (adapted from T Hengl)
setGeneric("buffer.dist", function(observations, predictionDomain, classes, width, ...) 
  standardGeneric("buffer.dist") )

setMethod("buffer.dist", signature(observations = "SpatialPointsDataFrame", predictionDomain = "SpatialPixelsDataFrame"), function(observations, predictionDomain, classes, width, ...){
  if(missing(width)){ width <- sqrt(areaSpatialGrid(predictionDomain)) }
  if(!length(classes)==length(observations)){ stop("Length of 'observations' and 'classes' does not match.") }
  ## remove classes without any points:
  xg = summary(classes, maxsum=length(levels(classes)))
  selg.levs = attr(xg, "names")[xg > 0]
  if(length(selg.levs)<length(levels(classes))){
    fclasses <- as.factor(classes)
    fclasses[which(!fclasses %in% selg.levs)] <- NA
    classes <- droplevels(fclasses)
  }
  ## derive buffer distances
  s <- list(NULL)
  for(i in 1:length(levels(classes))){
    s[[i]] <- raster::distance(rasterize(observations[which(classes==levels(classes)[i]),1]@coords, y=raster(predictionDomain)), width=width, ...)
  }
  s <- s[sapply(s, function(x){!is.null(x)})]
  s <- brick(s)
  s <- as(s, "SpatialPixelsDataFrame")
  s <- s[predictionDomain@grid.index,]
  return(s)
})

## function to approximate local homogeneity

library(datawizard)
loc_hom <- function(x){
  require(datawizard)
  y <- case_when(x>0 ~ 1/(x^2), x==0 ~ NA) 
  y1 <- case_when(is.na(y) ~ median(y, na.rm=TRUE), .default = y)
  y2 <- rescale(y1, to = c(1,6))
  return(y2)
}


## get data while not on MarVid
library(dplyr)
data <- read.csv(file.path(Data,"species_densities.csv"))
data <- data %>% rename(taxon=clean_taxonomy, density = density_n100m2) %>% dplyr::select(SampID, taxon, density)

sample_info <- read.csv(file.path(Data,"sample_info.csv"))
sample_info <- sample_info %>% rename(video.line=VL, section=SampID2, mean.latitude=y_coord, mean.longitude=x_coord) %>%
  dplyr::select(video.line, section, mean.latitude, mean.longitude, SampID)

