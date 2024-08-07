# ------------------------------------ #
# ---------- RARE REPTILES  ---------- #
# ------------------------------------ #

### Analysis fire effects on rarity and dominance in a reptile community 

# ---------- 01_process_data  ---------- #

### Script authors: Amber Lim & Annabel Smith 

# load workspace
load("04_workspaces/processed_data.RData")

# load packages
library("lme4"); library("vegan")

# -------- Pre-review data processing  -------- #

# ----
# Key questions:

# what can we tell about the influence of fire on rare species from comparing species richness to diversity estimates?

# how do diversity estimates which differ in their quantification of dominance and rarity compare? For example, is there stronger effects on the estimates which prioritise rarity?

# how does fire influence the overall abundance (pooled abundance) of 'rare' reptile species, classed as those with fewer than x captures?

# how does fire influence phylogenetic diversity? For example, are there a greater number of reptile families in unburnt mallee than burnt mallee? Does fire drive a phylogenetic simplification?

# how does fire influence functional rarity (sensu Violle et al. 2017, TREE)? We might not get this far... Just noting it down as an interesting aspect. 

# Following Rabinowitz (1981), we could look at rarity based on: geographic range, habitat specificity and local abundance. We've been talking about abundance-related rarity, but we could also look at geographical rarity - i.e. the proportion of sites in which each species occurs. Habitat specificity would have to be done through the literature. 

# THE DATA:

# full data set:

dat2<-read.table("01_data/all_reptile_data.txt",header=T)
head(dat2)

# species abundance data:
sp_div<-as.data.frame.matrix(t(table(dat2$name, dat2$site)))
head(sp_div,3); dim(sp_div)

# standardise the data by trap effort:
head(dat1, 3)
rownames(sp_div)

xdat<-dat1[which(dat1$site %in% rownames(sp_div)),]
xdat$site; xdat$season
unique(dat2$season)
head(xdat, 3); dim(xdat)

# combine seasons 5 and 6. Aggregate of effort by sites + seasons
eft<-data.frame(aggregate(xdat$effort~xdat$site, FUN=sum))
colnames(eft)<-c("site","effort")
eft

head(sp_div,3); dim(sp_div)
rownames(sp_div)==as.character(eft$site)

sdiv<-sp_div

# captures / 1000 trap nights. Scaling up to 2 seasons.
sp_div2<-round(apply(sdiv, 2, function(x) (x/eft$effort)*1000),0)
head(sp_div2, 3); dim(sp_div2)

# check that there are no zeros:
apply(sp_div2, 2, sum)

# site summaries:
sum_dat<-data.frame(site=rownames(sp_div2))

# total abundance:
sum_dat$abund<-rowSums(sp_div2)

# species richness:
sum_dat$sp_rich<-apply(sp_div2,1,function(x) length(which(x>0)))

head(sum_dat,3);dim(sum_dat)

# relative abundance:
rel_abund2<-apply(sp_div2, 2, function(x) x/sum_dat$abund)
head(rel_abund2,3); dim(rel_abund2)

# simpson's index:
sum_dat$simps_ind<-rowSums(rel_abund2^2)

# simpson's index:
sum_dat$simps_ind2<-1/rowSums(rel_abund2^2)

# evenness, from Simpson's:
sum_dat$even<-sum_dat$simps_ind2/sum_dat$sp_rich

# shannon's index:
log_relab<-log(rel_abund2)
log_relab<-t(apply(log_relab,1,function(x) ifelse(is.infinite(x)==T,0,x)))

sum_dat$shann_ind<--rowSums(rel_abund2*log_relab)

# evenness, from Shannon's:
sum_dat$Hmax<-log(sum_dat$sp_rich)
sum_dat$even2<-sum_dat$shann_ind/sum_dat$Hmax

# check simpson against vegan calc:
# plot(sum_dat$simps_ind,diversity(rel_abund2, index="simpson"))
cor.test(sum_dat$simps_ind2,diversity(rel_abund2, index="simpson"))

# check shannon against vegan calc:
# plot(sum_dat$shann_ind,diversity(rel_abund2, index="shannon"))
cor.test(sum_dat$shann_ind,diversity(rel_abund2, index="shannon"))

specnumber(rel_abund2)==sum_dat$sp_rich

sum_dat$fire_cat<-factor(c(rep("Unburnt",2),rep("Burnt",2),rep("Medium",2),rep("Burnt",2),rep("Unburnt",4),rep("Medium",2)),levels=c("Burnt","Medium","Unburnt"))
head(sum_dat); dim(sum_dat)

# species matrix into data frame for processing
sp_div2 <- as.data.frame(sp_div2)

# captures / 1000 trap nights. 2 seasons of data combined
head(sp_div2, 3); dim(sp_div2)

# rearranging so unburnt is baseline
sum_dat$fire_cat <- factor(sum_dat$fire_cat, levels = c("Unburnt", "Medium", "Burnt"))

# making site a factor
sum_dat$site <- as.factor(sum_dat$site)

# List of 14 sites, taking site column and calling "H" to create 2 different variables for site, either "H" or "P"
sum_dat$location <- NA
sum_dat$location[which(unlist(gregexpr("H", sum_dat$site)) == 1)] <- "H"
sum_dat$location[which(unlist(gregexpr("H", sum_dat$site)) == -1)] <- "P"

# processed data of 14 sites with species diversity metrics
head(sum_dat, 6);dim(sum_dat)

#list of species
colnames(sp_div2)

#summing columns
sum(sp_div2$A_nor)
sp_sum <- data.frame(sp=names(apply(sp_div2, MARGIN = 2, FUN = sum)),abundance=apply(sp_div2, MARGIN = 2, FUN = sum))
row.names(sp_sum) <- 1:nrow(sp_sum)
#write.table(sp_sum, file="SpeciesAbundance.txt", row.names = F, sep = "\t", quote = F)

head(sp_sum);dim(sp_sum)

# reading in species data, quantifying based on proportion of species (25%) vs maximum (5%)
dir()
sp_sum2 <- read.table("01_data/R_species_list.txt", header =T)

head(sp_sum2);dim(sp_sum2)

sum_dat$location<-NA
sum_dat$location[grep("H",sum_dat$site)]<-"Hincks"
sum_dat$location[grep("P",sum_dat$site)]<-"Pinks"
sum_dat$location <- as.factor(sum_dat$location)
sum_dat$fire_cat <- factor(sum_dat$fire_cat, levels = c("Unburnt", "Medium", "Burnt"))

head(sp_div2, 3); dim(sp_div2)

sp_25 <- sp_sum2$Species_Code[which(sp_sum2$sum_method_25 == "r")]
length(sp_25) # 14 rare species in the lowest 25% category

sp_5 <- sp_sum2$Species_Code[which(sp_sum2$max_5 == "r")]
length(sp_5) # 25 rare species in the max. 5% cateogry

# standardised abundance data for rare 25% species
sp_div_25 <- sp_div2[,which(colnames(sp_div2)%in%sp_25)]
head(sp_div_25, 3); dim(sp_div_25)

# standardised abundance data set for rare 5% species
sp_div_5 <- sp_div2[,which(colnames(sp_div2)%in%sp_5)]
head(sp_div_5, 3); dim(sp_div_5)

rare.data <- data.frame(site=names(rowSums(sp_div_25)), rel_abund_25 = rowSums(sp_div_25), pres_25 = ifelse(rowSums(sp_div_25) == 0, 0, 1), rel_abund_5 = rowSums(sp_div_5), pres_5 = ifelse(rowSums(sp_div_5) == 0, 0, 1))

row.names(rare.data) <- 1:nrow(rare.data)
head(sum_dat, 6);dim(sum_dat)
head(rare.data); dim(rare.data)

### UP TO HERE, NEED TO RE-DO THE SPECIES RICHNESS OF RARE SPECIES, BEFORE GOING TO BERGER PARKER

# What about species richness of rare species?

head(sp_div_5, 3); dim(sp_div_5)
head(sp_div_25, 3); dim(sp_div_25)

sr.rare<-data.frame(site=names(apply(sp_div_25,1,function(x) length(which(x>0)))),sr_25=apply(sp_div_25,1,function(x) length(which(x>0))),sr_5=apply(sp_div_5,1,function(x) length(which(x>0))))
rownames(sr.rare)<-1:nrow(sr.rare)

# Check that the sites are in the same order (these should all be TRUE):
table(sum_dat$site == sr.rare$site)
head(sr.rare); dim(sr.rare)

# combining data sets into sum_dat
sum_dat <- merge(sum_dat, rare.data, by="site", all.x = T, all.y = F)

sum_dat <- merge(sum_dat, sr.rare, by="site", all.x = T, all.y = F)
head(sum_dat, 6);dim(sum_dat)

# From the MS: "For the community of rare species, we did not model diversity metrics that required proportional abundances (Simpson’s Diversity, Shannon’s Diversity, Berger-Parker Dominance) because the abundance data were dominated by zeros and ones."

dir()
# save.image("04_workspaces/processed_data.RData")

# ----

# Peer-review processing updates

# R2 Comment 6 Line 363: There are no explicit models of dominance, such as the Berger-Parker index. Consider exploring Fisher’s alpha as an additional metric. 

# R2 Comment 13 A further way to explore how patterns of commonness and rarity reflect patterns of richness would be to rank the community from most to least common, and vice versa, and iteratively assess how well the richness calculated via these two ordered rankings matches patterns of overall species richness. This approach is similar to the method used by Lennon et al. https://onlinelibrary.wiley.com/doi/abs/10.1046/j.1461-0248.2004.00548.x

# -------- Berger-Parker  -------- #

# R2 Comment 6 Line 363: There are no explicit models of dominance, such as the Berger-Parker index. Consider exploring Fisher’s alpha as an additional metric. 

# Caruso et al. 2007: 
# "Equation (2) was thus used to calculate the Berger–Parker index for each replicate sample with no further standardization
# Eqn. 2:
# K = NS / NT
# where Ni is the abundance of the i-th species, NS indicates the abundance of the most abundant species, NT and S respectively are the total number of individuals and of species in a community sample..."

# Berger & Parker 1970"
# "...species dominance is given as:
# Dd = p max

# where p max is the maximum proportion of any one species in a sample..."

# Kitikidou et al. 2024:
# "The Berger–Parker index provides a simple quantification of the most abundant species in a given sample. It ranges from 0 to 1, with 0 indicating complete evenness or equal abundance among all species in the community. Conversely, a value of 1 indicates complete dominance, with one species accounting for all individuals in the community. 
# The Berger–Parker index is a simple and straightforward measure that can be easily calculated, making it useful for preliminary analyses or comparisons of dominance across different samples. However, it is important to note that the Berger–Parker index only considers the dominance of a single species and does not provide information about the overall diversity or richness of the community..."
"BP = max(pi)"
# The closer BP is to 0, the greater the diversity."

# It's the maximum proportional abundance


# Full data set (site x species matrix):
# 14 sites, standardised for trap effort (captures / 1000 trap nights)
head(sp_div2, 3); dim(sp_div2)

# Rare species(site x species matrix):
head(sp_div_25, 3); dim(sp_div_25)
head(sp_div_5, 3); dim(sp_div_5)

# Processed metrics with site data:
head(sum_dat, 2);dim(sum_dat)


# merge data:
sum_dat <- merge(sum_dat, sr.rare, by="site", all.x = T, all.y = F)
head(sum_dat, 6);dim(sum_dat)

# Relative abundance already calculated
head(rel_abund2,3); dim(rel_abund2)

bergpark<-data.frame(site=names(apply(rel_abund2,1,function(x) max(x))),bp_ind=apply(rel_abund2,1,function(x) max(x)))

bergpark$site==sum_dat$site

plot(bergpark$bp_ind, sum_dat$shann_ind)
plot(sum_dat$fire_cat, bergpark$bp_ind)



# Fisher's alpha:

fa.all<-fisher.alpha(sp_div2)
fa.sp25<-fisher.alpha(sp_div_25)
fa.sp5<-fisher.alpha(sp_div_5)

