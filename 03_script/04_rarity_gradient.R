# ------------------------------------ #
# ---------- RARE REPTILES  ---------- #
# ------------------------------------ #

### Analysis fire effects on rarity and dominance in a reptile community 

# ---------- 04_rarity_gradient  ---------- #

### Script authors: Annabel Smith & Amber Lim

# If starting fresh:
# load("04_workspaces/analysed_data.RData"); rm(list = setdiff(ls(), c("sum_dat", "dat1","dat2", "sp_div","sp_div2","rel_abund2","sp_sum2")))

# Load workspace: 
load("04_workspaces/rarity_gradient.RData")

# save.image("04_workspaces/rarity_gradient.RData")

# Load functions:
invisible(lapply(paste("02_functions/",dir("02_functions"),sep=""), function(x) source(x)))

# THE DATA:

# diversity data fully processed and ready to analyse:
head(sum_dat, 3);dim(sum_dat)

# dat1 from the original analysis; raw, unstandardised abundances with seasons separate:
head(dat1,3); dim(dat1)

# full, unprocessed individual reptile records:
head(dat2); dim(dat2)

# species abundance data (unstandardised):
head(sp_div,3); dim(sp_div)

# relative abundance:
head(rel_abund2,3); dim(rel_abund2)

# captures / 1000 trap nights:
head(sp_div2, 3); dim(sp_div2)

# species data 
head(sp_sum2);dim(sp_sum2)

# Check that the abundances are correct:
sp_sum_species<-sp_sum2$Abundance
names(sp_sum_species)<-sp_sum2$Species_Code

# These should all be TRUE (41 species)
table(sp_sum_species[order(sp_sum_species)]==colSums(sp_div2)[order(colSums(sp_div2))])

# Abundance data:
ab_dat<-sp_sum2[,c("Species_Code","Abundance")]
ab_dat<-ab_dat[order(ab_dat$Abundance),]
rownames(ab_dat)<-1:nrow(ab_dat)

# Abundance ranked from low to high
head(ab_dat, 3);dim(ab_dat)

# Species richness by site:
sr_dat<-sum_dat[,c("site","abund","sp_rich","fire_cat","location")]

# Standardised abundance for calculating species richness:
head(sp_div2, 3); dim(sp_div2)

# Check species richness for all species:
rich.all<-apply(sp_div2,1,function(x) length(which(x>0)))

# These should all be TRUE (14 sites):
table(sr_dat$sp_rich==rich.all)


# Rank species from low to high

head(ab_dat, 3);dim(ab_dat) # species total abundances
head(sr_dat, 3);dim(sr_dat) # summarised species richness
head(sp_div2, 3); dim(sp_div2) # species abundances by site

low_to_high<-list()

for (i in 1:length(ab_dat$Species_Code)){
  
  sp.thisrun<-ab_dat$Species_Code[i]
  all.sp.thisrun<-ab_dat$Species_Code[1:i]
  dat.thisrun<-data.frame(sp_div2[all.sp.thisrun])
  rich.thisrun<-data.frame(apply(dat.thisrun,1,function(x) length(which(x>0))))
  colnames(rich.thisrun)<-paste("lh",i,sep="")
  low_to_high[[i]]<-rich.thisrun
  
} # close for species

low_high <- data.frame(do.call(cbind, low_to_high))
low_high$site<-rownames(low_high)
head(low_high[,1:12], 3); dim(low_high)
head(low_high[,(ncol(low_high)-10):ncol(low_high)], 3); dim(low_high)

# Add ranked richness to main data frame:

sr_dat<-merge(sr_dat, low_high, by="site", all.x=T, all.y=F)

# The original species richness should be the same as the highest calculation:
table(sr_dat$sp_rich==sr_dat$lh41)

head(sr_dat, 3);dim(sr_dat)

dev.new(width=8,height=12,dpi=70,pointsize=16, noRStudioGD = T)
par(mfrow=c(7,6), mgp=c(3,1,0), mar=c(4,2,1,0))

lhs<-colnames(sr_dat)[grep("lh",colnames(sr_dat))]

for(i in 1:length(lhs)){
  
  lh.thisrun<-lhs[i]
  ind.thisrun<-which(colnames(sr_dat)==lh.thisrun)
  plot(sr_dat$fire_cat, sr_dat[,ind.thisrun], las=1, ylab="", xlab="")

} # close i lhs
  
  
# save.image("04_workspaces/rarity_gradient.RData")







