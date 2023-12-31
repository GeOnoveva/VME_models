library(dplyr)
library(caret)

vme <- "Hard bottom sponge aggregations"
sp = taxonary %>% filter(VME_Burgos_etal_2020==vme) %>% select(Reference_List) %>% pull
e <- read.csv(file.path(savepath, "e.csv"))

# obtain only rows of interest
data_extract <- data %>% filter(taxon%in%sp) %>% left_join(e) # unexpected many-to-many relationships found!!!
  
# create empty threshold vector
threshold <- vector(mode = "numeric",
                   length = length(sp))
  
for (i in 1:length(sp)){
    
  threshold[i] <- (data %>% filter(taxon==sp[i]) %>% pull(density) %>% sort %>%  quantile(.))[3]
    
}
  
thresholdt <- data.frame(cbind(sp, threshold))
 
# for each row check that the density is greater than threshold. 
# Result should be vector of TRUEs and FALSEs
  
data_extract_merged <-merge(data_extract,thresholdt,by.x="taxon",by.y="sp",all.x=TRUE)
data_extract_merged_reduced <- data_extract_merged %>% filter(density>threshold) %>%
  select(taxon,          
         bathy,                        
         BO22_carbonphytoltmax_bdmean,
         BO22_carbonphytoltmax_ss,
         BO22_carbonphytoltmin_bdmean,
         BO22_carbonphytoltmin_ss,  
         BO22_carbonphytomean_bdmean,
         BO22_carbonphytomean_ss,
         BO22_carbonphytorange_bdmean,
         BO22_carbonphytorange_ss,
         BO22_chloltmax_bdmean,
         BO22_chloltmax_ss,    
         BO22_chloltmin_bdmean,       
         BO22_chloltmin_ss,  
         BO22_chlomean_bdmean,    
         BO22_chlomean_ss,  
         BO22_chlorange_bdmean,       
         BO22_chlorange_ss,        
         BO22_dissoxltmax_bdmean,   
         BO22_dissoxltmin_bdmean,  
         BO22_dissoxmean_bdmean,      
         BO22_dissoxrange_bdmean,   
         BO22_icecoverltmax_ss,     
         BO22_icecoverltmin_ss,     
         BO22_icecovermean_ss,        
         BO22_icecoverrange_ss,     
         BO22_icethickltmax_ss,     
         BO22_icethickltmin_ss,    
         BO22_icethickmean_ss,        
         BO22_icethickrange_ss,    
         BO22_ironltmax_bdmean,     
         BO22_ironltmin_bdmean,     
         BO22_ironmean_bdmean,        
         BO22_ironrange_bdmean,     
         BO22_lightbotltmax_bdmean,  
         BO22_lightbotltmin_bdmean,  
         BO22_lightbotmean_bdmean,    
         BO22_lightbotrange_bdmean,  
         BO22_nitrateltmax_bdmean,   
         BO22_nitrateltmin_bdmean, 
         BO22_nitratemean_bdmean,     
         BO22_nitraterange_bdmean,    
         BO22_phosphateltmax_bdmean, 
         BO22_phosphateltmin_bdmean,  
         BO22_phosphatemean_bdmean,   
         BO22_phosphaterange_bdmean, 
         BO22_ppltmax_bdmean,       
         BO22_ppltmax_ss,            
         BO22_ppltmin_bdmean,         
         BO22_ppltmin_ss,            
         BO22_ppmean_bdmean,        
         BO22_ppmean_ss,           
         BO22_pprange_bdmean,         
         BO22_pprange_ss,           
         BO22_silicateltmax_bdmean,   
         BO22_silicateltmin_bdmean,   
         BO22_silicatemean_bdmean,    
         BO22_silicaterange_bdmean,   
         CDirmax_Robinson,         
         CDirmean_Robinson,        
         CDirmin_Robinson,            
         CDirsd_Robinson,          
         cobB,                    
         CSpdmax_Robinson,        
         CSpdmean_Robinson,           
         CSpdmin_Robinson,         
         CSpdsd_Robinson,          
         diffME3,                   
         diffME9,                     
         gravel,                      
         landscape,                  
         MLDmax_Robinson,             
         MLDmean_Robinson,            
         MLDmin_Robinson,             
         MLDsd_Robinson,             
         msr1_mag,                    
         msr5_mag,                    
         MS_biogeo05_dist_shore_5m,   
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
         Smax_Robinson,              
         Smean_Robinson,             
         Smin_Robinson,               
         spd_max,                     
         spd_mean,                   
         spd_min,                  
         spd_std,                     
         Ssd_Robinson,              
         temp_max,                 
         temp_mean,                
         temp_min,                    
         temp_std,               
         Tmax_Robinson,           
         Tmean_Robinson,       
         Tmin_Robinson,               
         Tsd_Robinson,   
         Umax_Robinson,    
         Umean_Robinson,   
         Umin_Robinson,               
         Usd_Robinson,   
         u_bott_mean,     
         Vmax_Robinson,    
         Vmean_Robinson,              
         Vmin_Robinson,    
         Vsd_Robinson,      
         v_bott_mean    )

data_extract_merged_reduced <- data_extract_merged_reduced %>% filter(complete.cases(.))

library(corrplot)

df <- e %>% dplyr::select(-c(1,119:123)) %>% dplyr::select(!BO22_icecoverltmin_ss) %>% filter(complete.cases(.))


# remove zero variance predictors
zv <- apply(df, 2, function(x) length(unique(x)) == 1)
dfr <- df[, !zv]
n=length(colnames(dfr))

# remove highly colinear

correlationMatrix <- cor(dfr[,1:n],use="complete.obs")
namesToDrop <- findCorrelation(cor(dfr, use="pairwise.complete.obs"), cutoff = 0.8, names=TRUE)

data_extract_merged_reduced_new <- data_extract_merged_reduced %>%
  select(-na.omit(match(namesToDrop, colnames(data_extract_merged_reduced))))

