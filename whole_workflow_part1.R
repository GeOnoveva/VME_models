# Requirements
## Taxonary
## data
## sample info
## tab station
## e

# Set up parallel processing
#cl <- makeCluster(39, type='PSOCK') # leave 25 cores available to other people!
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)

# Set vme
vme <- "Deep-sea sea pen communities"
remove <- "nothing"
indicators <- taxonary$Reference_List[which(taxonary$VME_Burgos_etal_2020==vme)]
#indicators <- indicators[-which(indicators%in%remove)]

# Filter data by VME indicators
resp1 <- data %>% group_by(SampID) %>%
  summarize(tot_rich = length(SampID), tot_dens = sum(density))
resp2 <- data %>% filter(taxon%in%indicators) %>%
  group_by(SampID) %>%
  summarize(vmeind_density = sum(density))
resp <- left_join(resp1,resp2) %>% mutate(vmeind_density=replace_na(vmeind_density,0))

# Save output
rm(resp1, resp2)
write.csv(resp, file.path(savepath, paste(vme, "minus", paste(substr(remove, 1,7), collapse=", "), "resp.csv", collapse=" ")))

# Aggregate to whole video line and calculate sd
respw <- resp %>% mutate(video.line=sub("_.*", "", SampID))
respw <- respw %>% group_by(video.line) %>%
  summarize(mean_rich = mean(tot_rich), mean_dens = mean(vmeind_density), sd_dens = sd(vmeind_density)) %>%
  mutate(video.line = as.integer(video.line))

# Make spatial object
resp_spat1 <- resp %>% mutate(mean_long=as.numeric(pull(
  dplyr::filter(data.frame(lookup(SampID,dplyr::select(sample_info, c(SampID, mean.longitude)), missing = NULL)), 
                !duplicated(data.frame(qdapTools::lookup(SampID,dplyr::select(sample_info, c(SampID, mean.longitude)), missing = NULL)))
  ))),
  mean_lat=as.numeric(pull(
    dplyr::filter(data.frame(lookup(SampID,dplyr::select(sample_info, c(SampID, mean.latitude)), missing = NULL)), 
                  !duplicated(data.frame(qdapTools::lookup(SampID,dplyr::select(sample_info, c(SampID, mean.latitude)), missing = NULL)))
    )))) %>%
  filter(!is.na(mean_long))
resp_spat <- SpatialPointsDataFrame(dplyr::select(resp_spat1, c(mean_long, mean_lat)), dplyr::select(resp_spat1, c(SampID,tot_rich, tot_dens, vmeind_density)))  
#%>% utmize
writeOGR(resp_spat, savepath, paste(vme, "minus", paste(substr(remove, 1,7), collapse=", "), "resp_spat"), driver = "ESRI Shapefile", overwrite_layer = TRUE)

# make resp have the same number of rows as resp_spat:
resp <- resp %>% filter(SampID%in%resp_spat@data$SampID)

# whole video line level
resp_spat1 <- respw %>% mutate(mean_long=as.numeric(pull(
  dplyr::filter(data.frame(lookup(video.line,dplyr::select(tab_station, c(sample_no, lon_mid_dec)), missing = NULL))))),
  mean_lat=as.numeric(pull(
    dplyr::filter(data.frame(lookup(video.line,dplyr::select(tab_station, c(sample_no, lat_mid_dec)), missing = NULL)))))) %>%
  dplyr::filter(!is.na(mean_long)) %>% rename(vmeind_density = mean_dens)
coords <- coordinates(dplyr::select(resp_spat1, c(mean_long, mean_lat)))
resp_spat_vl <- SpatialPointsDataFrame(coords, dplyr::select(resp_spat1, c(video.line,mean_rich, vmeind_density, sd_dens))) %>%
  utmize
writeOGR(resp_spat_vl, savepath, paste(vme, "minus", paste(substr(remove, 1,7), collapse=", "),  "resp_spat_vl"), driver = "ESRI Shapefile", overwrite_layer = TRUE)

# make resp have the same number of rows as resp_spat:
respw <- respw %>% filter(video.line%in%resp_spat_vl@data$video.line)

# Save outputs
rm(resp_spat1)
write.csv(resp, file.path(savepath, paste(vme,"minus", paste(substr(remove, 1,7), collapse=", "),"resp new.csv", collapse=" ")))
write.csv(respw, file.path(savepath, paste(vme,"minus", paste(substr(remove, 1,7), collapse=", "), "st_dev_vme_dens.csv")), row.names = FALSE)
#saveRDS(resp_spat_vl, file.path(savepath, paste(vme,"minus", paste(substr(remove, 1,7), collapse=", "), "resp_Spat_vl.rds")))
saveRDS(resp_spat, file.path(savepath, paste(vme,"minus", paste(substr(remove, 1,7), collapse=", "), "resp_spat.rds")))
stdev <- respw

# Inputs
v <- left_join(resp,data.frame(e), by="SampID")

# join standard deviation of vme density
v <- v %>% left_join(stdev, by="video.line")

# save
write.csv(dplyr::select(v,
                        SampID,
                        tot_rich,
                        tot_dens,
                        vmeind_density,
                        mean_dens,
                        mean_rich,
                        sd_dens,
                        bathy,
                        BO22_carbonphytoltmax_bdmean,
                        BO22_carbonphytoltmax_ss,
                        BO22_carbonphytoltmin_bdmean,
                        BO22_carbonphytoltmin_ss,
                        BO22_carbonphytomean_bdmean,
                        BO22_carbonphytomean_ss,
                        BO22_carbonphytorange_bdmean,
                        BO22_carbonphytorange_ss,    
                        BO22_chloltmin_bdmean,    
                        BO22_chlomean_ss,    
                        BO22_dissoxltmax_bdmean,    
                        BO22_dissoxrange_bdmean,    
                        BO22_icecovermean_ss,    
                        BO22_icethickltmin_ss,    
                        BO22_ironltmax_bdmean,    
                        BO22_ironrange_bdmean,    
                        BO22_lightbotmean_bdmean,    
                        BO22_nitrateltmin_bdmean,    
                        BO22_phosphateltmax_bdmean,    
                        BO22_phosphaterange_bdmean,    
                        BO22_ppltmin_bdmean,    
                        BO22_ppmean_ss,    
                        BO22_silicateltmax_bdmean,    
                        BO22_silicaterange_bdmean,    
                        BO22_chloltmax_bdmean,
                        BO22_chloltmax_ss,
                        BO22_chloltmin_ss,
                        BO22_chlomean_bdmean,
                        BO22_chlorange_bdmean,
                        BO22_chlorange_ss,
                        BO22_dissoxltmin_bdmean,
                        BO22_dissoxmean_bdmean,
                        BO22_icecoverltmax_ss,
                        BO22_icecoverltmin_ss,
                        BO22_icecoverrange_ss,
                        BO22_icethickltmax_ss,
                        BO22_icethickmean_ss,
                        BO22_icethickrange_ss,
                        BO22_ironltmin_bdmean,
                        BO22_ironmean_bdmean,
                        BO22_lightbotltmax_bdmean,
                        BO22_lightbotltmin_bdmean,
                        BO22_lightbotrange_bdmean,
                        BO22_nitrateltmax_bdmean,
                        BO22_nitratemean_bdmean,
                        BO22_nitraterange_bdmean,
                        BO22_phosphateltmin_bdmean,
                        BO22_phosphatemean_bdmean,
                        BO22_ppltmax_bdmean,
                        BO22_ppltmax_ss,
                        BO22_ppltmin_ss,
                        BO22_ppmean_bdmean,
                        BO22_pprange_bdmean,
                        BO22_pprange_ss,
                        BO22_silicateltmin_bdmean,
                        BO22_silicatemean_bdmean,
                        CDirmax_Robinson,
                        CDirmean_Robinson,
                        CDirmin_Robinson,    
                        CDirsd_Robinson,
                        CSpdmax_Robinson,    
                        CSpdmean_Robinson,
                        CSpdmin_Robinson,
                        CSpdsd_Robinson,    
                        MLDmax_Robinson,
                        MLDmean_Robinson,    
                        MLDmin_Robinson,
                        MLDsd_Robinson,
                        MS_biogeo05_dist_shore_5m,
                        Smax_Robinson,    
                        Smean_Robinson,
                        Smin_Robinson,
                        Ssd_Robinson,
                        Tmax_Robinson,                  
                        Tmean_Robinson,
                        Tmin_Robinson,
                        Tsd_Robinson,                  
                        Umax_Robinson,
                        Umean_Robinson,
                        Umin_Robinson,                  
                        Usd_Robinson,
                        Vmax_Robinson,                  
                        Vmean_Robinson,
                        Vmin_Robinson,
                        Vsd_Robinson,
                        cobB,
                        diffME3,
                        diffME9,
                        gravel,    
                        landscape,
                        msr1_mag,    
                        msr5_mag,
                        mud,    
                        rock,
                        salt_max,
                        salt_mean,    
                        salt_min,
                        salt_std,
                        sand,    
                        sedclass,
                        seddan,
                        sedmil,    
                        slope3,
                        slope9,
                        spd_max,                  
                        spd_mean,
                        spd_min,
                        spd_std,                  
                        temp_max,
                        temp_mean,                  
                        temp_min,
                        temp_std,
                        u_bott_mean,
                        v_bott_mean		    
), file.path(savepath, paste(vme, "minus", paste(substr(remove, 1,7), collapse=", "), "v.csv", collapse=" ")), row.names = FALSE)

v <- read.csv(file.path(savepath, paste(vme, "minus", paste(substr(remove, 1,7), collapse=", "), "v.csv", collapse=" ")))
v <- v %>% 
  mutate_at(vars(  landscape, sedclass, seddan, sedmil    ), factor)
v <- v %>% filter(complete.cases(.))

# Select environmental variables
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
a=which(colnames(v)=="bathy")
b=which(colnames(v)=="v_bott_mean")
c=which(colnames(v)=="vmeind_density")

# Configure multicore
registerDoParallel(cl)

# variable selection
results <- rfe(v[,a:b], v[,c], sizes=c(1:(b-a)), rfeControl=control)

# chosen features
selvar <- predictors(results)
write.csv(selvar, file.path(savepath, paste(vme, "minus", paste(substr(remove, 1,7), collapse=", "), "selvar.csv", collapse=" ")),
          row.names = FALSE)
selvar <- read.csv(file.path(savepath, paste(vme, "minus", paste(substr(remove, 1,7), collapse=", "), "selvar.csv")))
selvar <- pull(selvar)

all_names <- paste(#c(names(resp.dist0),
  selvar, collapse="+")
fmla0 <- as.formula(paste("vmeind_density ~ ", 
                          all_names))

# fine tune mtry parameter
control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
set.seed(1)
tunegrid <- expand.grid(.mtry=c(1:length(selvar)))
rf_gridsearch <- train(fmla0, data=v, method="qrf", metric="Rsquared", tuneGrid=tunegrid, trControl=control)
print(rf_gridsearch)
mtry=rf_gridsearch$bestTune
write.csv(mtry, file.path(savepath, paste(vme, "minus", paste(substr(remove, 1,7), collapse=", "), "mtry.csv", collapse=" ")),
          row.names = FALSE)

