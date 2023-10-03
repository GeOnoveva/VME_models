library(nicheROVER)
data(fish) # 4 fish, 3 isotopes
#aggregate(fish[2:4], fish[1], mean) # isotope means calculated for each species

nsamples <- 1000

#fish.par <- tapply(1:nrow(fish), fish$species,
#                   function(ii) niw.post(nsamples = nsamples, X = fish[ii,2:4]))


d<-data_extract_merged_reduced_new


zv <- apply(d, 2, function(x) length(unique(x)) == 1)
d <- d[, !zv]
d <- d %>% filter(!taxon%in%c( "Phakellia ventilabrum",
                                     "Polymastia grimaldii",
                                     "Tetilla sp."))

d <- d %>% select(-c(  landscape, sedclass, seddan, sedmil, BO22_icecovermean_ss, BO22_lightbotltmax_bdmean,
                             BO22_ppltmin_bdmean, BO22_ppltmin_ss))

d$taxon<- as.factor(as.character(d$taxon))


p <- ggplot(d, aes(x=sand         , color = taxon)) + 
  geom_density() +scale_color_brewer(palette="Set3")
p



d.par <- tapply(1:nrow(d), d$taxon,
                   function(ii) niw.post(nsamples = nsamples, X = d[ii,2: dim(d)[2]]))

# Overlap calculation.  use nsamples = nprob = 10000 (1e4) for higher accuracy.
# the variable over.stat can be supplied directly to the overlap.plot function

over.stat <- overlap(d.par, nreps = nsamples, nprob = 1e3, alpha = c(.95, 0.99))

#The mean overlap metrics calculated across iteratations for both niche 
#region sizes (alpha = .95 and alpha = .99) can be calculated and displayed in an array.
over.mean <- apply(over.stat, c(1:2,4), mean)*100

coex<-round(over.mean, 2)[,,1] # for alpha=95%

result <- cbind(rowMeans(data.frame(coex), na.rm=TRUE), colMeans(data.frame(coex), na.rm=TRUE)) %>%
  data.frame 

rogue <- filter(result, X1<10 & X2<10) %>% row.names(.)
rogue
