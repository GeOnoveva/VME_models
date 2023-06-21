
# Remember to wipe the environment when running from the top!

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc()

# Paths:

## Path to local downloads folder, where MarVid reports get saved
Data <- "C:\\Users\\genoveva\\Downloads"

## Path to HI OneDrive folder, where Taxonary resides
SharedData <- "C:\\Users\\genoveva\\Havforskningsinstituttet\\Video analysis lab - Reference List\\"

## Path to ftp site, where environmental layers (already cropped and aligned) reside
ftp_url <- "ftp://ftp2.ngu.no/toGeno/EnvLayersNiN2022/all"
ftp_filepath <- paste0(strsplit(ftp_url, "\\//")[[1]][1], "//havforsk:Havf8776@", 
                   strsplit(ftp_url, "\\//")[[1]][2])#insert the credentials into the url

## Path to folder where to save intermediate objects created along the way and needed in subsequent steps (you need to create this folder if it doesn't exist)
savepath <- "..\\..\\Intermediate_objects"

## Final products to be shared:
outpath <- "U:\\Mareano\\Sårbare_habitater_kart"

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

### function to remove from v observations where response in the low class, then split the rest randomly, all while keeping lines together

stratisplit <- function (x, n, keep, p, video.line) {
  x<-v
  n<-3
  keep=0.25
  p=0.3
  video.line=video.line # lookup video line by matching v to sample info on SampID
  require(BAMMtools, groupdata2, tidyverse)
  
  # stratification by response
  densclas <- cut(v$vmeind_density, breaks = getJenksBreaks(v$vmeind_density, n+1),# 3 strata (jenks breaks)
                  include.lowest = TRUE,
                  labels = FALSE) 
  
  # first, find out which video lines have most sections in the low dens class: defined as having max 1 section in a higher class:
  d1 <- data.frame(cbind(densclas,video.line))
  d1$value <- 1
  d2 <- pivot_wider(d1, names_from = "densclas", values_fn = sum, values_fill = 0)
  colnames(d2)[2:4] <- c("low", "interm", "high")
  d2 <- d2 %>% mutate(d2, tot = rowSums(select(.,low, interm, high)),
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
    select(-which(colnames(.)==".partitions")) %>% rbind(.,d4)
  return(list(v_test, v_train))
}

## function to polygonize the high density predicted areas

threspol <- function(thresvalue, raster, min.area, max.hole){
  rcl <- matrix(data = c(0,thresvalue,thresvalue,max(values(raster), na.rm =TRUE),0,1), nrow=2)
  rec <- reclassify(raster, rcl, include.lowest = TRUE)
  pol <- rasterToPolygons(rec, fun=function(x){x==1}, dissolve = TRUE)
  pol <- as(spatialEco::explode(pol), "Spatial")
  pol$area_sqkm <- area(pol) / 1000000
  pol <- pol[pol$area_sqkm>min.area,]
  area_thresh <- units::set_units(max.hole, km^2)
  pol <- fill_holes(pol, threshold = area_thresh)
  pol <- smoothr::smooth(pol, method = "ksmooth")
  pol_dens <- pol
  return(pol)
}
