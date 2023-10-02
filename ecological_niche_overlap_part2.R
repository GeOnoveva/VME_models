library(nicheROVER)
data(fish) # 4 fish, 3 isotopes
#aggregate(fish[2:4], fish[1], mean) # isotope means calculated for each species
nsamples <- 1000

fish.par <- tapply(1:nrow(fish), fish$species,
                   function(ii) niw.post(nsamples = nsamples, X = fish[ii,2:4]))


fish<-data_extract_merged_reduced
fish$taxon<- as.factor(as.character(fish$taxon))

# remove colinear variables

fish <- fish %>% filter(!taxon%in%c("Polymastia grimaldii", "Phakellia ventilabrum",
                                    "Tetilla sp.", "Craniella sp."))

fish.par <- tapply(1:nrow(fish), fish$taxon,
                   function(ii) niw.post(nsamples = nsamples, X = fish[ii,2:116]))

# Overlap calculation.  use nsamples = nprob = 10000 (1e4) for higher accuracy.
# the variable over.stat can be supplied directly to the overlap.plot function

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1e3, alpha = c(.95, 0.99))

#The mean overlap metrics calculated across iteratations for both niche 
#region sizes (alpha = .95 and alpha = .99) can be calculated and displayed in an array.
over.mean <- apply(over.stat, c(1:2,4), mean)*100
round(over.mean, 2)




# load pred out environment, which contains e

