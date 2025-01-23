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

# Load libraries:
library("MASS")
library(AICcmodavg)

# THE DATA ----

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

# ----

head(ab_dat, 3);dim(ab_dat) # species total abundances
head(sr_dat, 3);dim(sr_dat) # summarised species richness
head(sp_div2, 3); dim(sp_div2) # species abundances by site

# Rank species along rarity gradients ----

# low to high

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

sr_lh<-merge(sr_dat, low_high, by="site", all.x=T, all.y=F)

# The original species richness should be the same as the highest calculation:
table(sr_dat$sp_rich==sr_lh$lh41)

# high to low

ab_dat2<-ab_dat[order(ab_dat$Abundance, decreasing = T),]
rownames(ab_dat2)<-1:nrow(ab_dat2)
head(ab_dat2, 3);dim(ab_dat2)

high_to_low<-list()

for (i in 1:length(ab_dat2$Species_Code)){
  
  sp.thisrun<-ab_dat2$Species_Code[i]
  all.sp.thisrun<-ab_dat2$Species_Code[1:i]
  dat.thisrun<-data.frame(sp_div2[all.sp.thisrun])
  rich.thisrun<-data.frame(apply(dat.thisrun,1,function(x) length(which(x>0))))
  colnames(rich.thisrun)<-paste("hl",i,sep="")
  high_to_low[[i]]<-rich.thisrun
  
} # close for species

high_low <- data.frame(do.call(cbind, high_to_low))
high_low$site<-rownames(high_low)
head(high_low[,1:12], 3); dim(high_low)
head(high_low[,(ncol(high_low)-10):ncol(high_low)], 3); dim(high_low)

# Add ranked richness to main data frame:

sr_hl<-merge(sr_dat, high_low, by="site", all.x=T, all.y=F)

# The original species richness should be the same as the highest calculation:
table(sr_dat$sp_rich==sr_hl$hl41)

# ----

# Model fire effect along rarity gradient ----

# low to high:
lhs<-colnames(sr_lh)[grep("lh",colnames(sr_lh))]

# model fire effect for each lh value
lhs_dat<-sr_lh[,lhs]
colSums(lhs_dat)

head(sr_lh[,1:15], 3);dim(sr_lh)

lh_mods<-list()

lh_res<-data.frame(lh=lhs,total_sr=colSums(lhs_dat),AICc_fire=NA, AICc_null=NA)
rownames(lh_res)<-1:nrow(lh_res)

for (i in 1:length(lhs)){
  
  resp.thisrun<-lhs[i]
  m1_fire <- glm.nb(get(resp.thisrun)~fire_cat,data = sr_lh)
  m1_null <- glm.nb(get(resp.thisrun)~1,data = sr_lh)
  summary(m1_fire)
  
  lh_res$AICc_fire[i]<-AICc(m1_fire)
  lh_res$AICc_null[i]<-AICc(m1_null)
  
  lh_mods[[i]]<-summary(m1_fire)$coefficients
  
} # close run models

# Calculate Delta AICc
lh_res$delta<-lh_res$AICc_fire-lh_res$AICc_null
head(lh_res,3); dim(lh_res)


# high to low:
hls<-colnames(sr_hl)[grep("hl",colnames(sr_hl))]

# model fire effect for each lh value
hls_dat<-sr_hl[,hls]
colSums(hls_dat)

head(sr_hl[,1:15], 3);dim(sr_hl)

hl_mods<-list()

hl_res<-data.frame(hl=hls,total_sr=colSums(hls_dat),AICc_fire=NA, AICc_null=NA)
rownames(hl_res)<-1:nrow(hl_res)

for (i in 1:length(hls)){
  
  resp.thisrun<-hls[i]
  m1_fire <- glm.nb(get(resp.thisrun)~fire_cat,data = sr_hl)
  m1_null <- glm.nb(get(resp.thisrun)~1,data = sr_hl)
  summary(m1_fire)
  
  hl_res$AICc_fire[i]<-AICc(m1_fire)
  hl_res$AICc_null[i]<-AICc(m1_null)
  
  hl_mods[[i]]<-summary(m1_fire)$coefficients
  
} # close run models

# Calculate Delta AICc
hl_res$delta<-hl_res$AICc_fire-hl_res$AICc_null
head(hl_res,3); dim(hl_res)

# ----

# Plot delta AICc ----

dev.new(width=6,height=4,dpi=70,pointsize=16, noRStudioGD = T)
par(mfrow=c(1,1), mgp=c(2.2,0.8,0), mar=c(3.5,3.5,1,1), oma=c(0,0,0,7))
plot(1:nrow(lh_res), lh_res$delta, type="l", xlab="Rarity gradient", ylab="", las=1)
title(ylab=as.expression(bquote(Delta~"AICc")))
lines(0:45,rep(2,46), col="red",lty=1, lwd=1.2)

lines(1:nrow(hl_res), hl_res$delta, col="cornflowerblue",lwd=1.2)

# plot legend:
par(xpd=NA)
legend(x=44,y=8, title = "Species added", legend = c("Low to high", "High to low"), lwd = 1.2, col = c("black","cornflowerblue"), bty = "n", title.adj=0)
par(xpd=F)

# ----





# Box plots of raw data:
dev.new(width=8,height=12,dpi=70,pointsize=16, noRStudioGD = T)
par(mfrow=c(7,6), mgp=c(3,1,0), mar=c(4,2,1,0))

for(i in 1:length(lhs)){
  
  lh.thisrun<-lhs[i]
  ind.thisrun<-which(colnames(sr_dat)==lh.thisrun)
  plot(sr_dat$fire_cat, sr_dat[,ind.thisrun], las=1, ylab="", xlab="")

} # close i box plot
  
  
# save.image("04_workspaces/rarity_gradient.RData")







