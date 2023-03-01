# Ideas for Amber's project
# Feb 2023

# what can we tell about the influence of fire on rare species from comparing species richness to diversity estimates?

# how do diversity estimates which differ in their quantification of dominance and rarity compare? For example, is there stronger effects on the estimates which prioritise rarity?

# how does fire influence the overall abundance (pooled abundance) of 'rare' reptile species, classed as those with fewer than x captures?

# how does fire influence phylogenetic diversity? For example, are there a greater number of reptile families in unburnt mallee than burnt mallee? Does fire drive a phylogenetic simplification?

# how does fire influence functional rarity (sensu Violle et al. 2017, TREE)? We might not get this far... Just noting it down as an interesting aspect. 

# Following Rabinowitz (1981), we could look at rarity based on: geographic range, habitat specificity and local abundance. We've been talking about abundance-related rarity, but we could also look at geographical rarity - i.e. the proportion of sites in which each species occurs. Habitat specificity would have to be done through the literature. 

#testing amber 1 march


# For community ecology lecture
# March 2020

# Set up the data:

dat1<-read.table("sp_summaries.txt",header=T)
dat1$TSF2<-dat1$TSF^2
dat1$TSFscale<-scale(dat1$TSF,center=T,scale=T)
dat1$TSF2scale<-scale(dat1$TSF^2,center=T,scale=T)
dat1$obs_effect<-as.numeric(1:nrow(dat1))
dat1$effort<-as.numeric(dat1$effort)
head(dat1); dim(dat1)

# How many captures per location for each species?

spnames<-names(dat1[,9:(length(dat1)-4)])
names_lgth<-length(spnames)
caps<-matrix(data=NA,nrow=length(spnames),ncol=2)

for(i in 1:names_lgth){
  
  cols_thisrun<-colnames(dat1)==spnames[i]
  caps[i,]<-tapply(dat1[,cols_thisrun],dat1$location,sum, na.rm=T)
}
caps<-data.frame(sp=spnames,caps)
colnames(caps)<-c("sp",levels(dat1$location))
caps


# Set up a factor level which assigns each species to one of three groups: 1 = Hincks only, 2 = Pinks only, and 3 = both. The cut-off for analysis is >=28 Hincks and >=38 Pinks.

sites_tab<-data.frame(Hsites=length(which(table(dat1$site, dat1$location)[,1]>0)),Psites=length(which(table(dat1$site, dat1$location)[,2]>0)))

spnames2<-caps[,1]
lgth2<-length(caps[,1])
new.x<-matrix(data=NA,nrow=lgth2,ncol=1)

for (i in 1:lgth2){
  
  new.x[i,]<- if (caps[i,2]>=28 & caps[i,3]<38) 1 
  else if (caps[i,2]<28 & caps[i,3]>=38) 2
  else 3
}
caps$mod_type<-as.factor(new.x)
caps



# STANDARDISE THE DATA BY TRAP EFFORT (dat4). 

# This will turn the data into something like "the number of captures per 500 trap nights".

head(dat1)
lgth2<-length(caps[,1])
lgth5<-length(dat1$obs)
sp_cols<-dat1[,9:(length(dat1)-4)]
div<-matrix(data=NA, nrow=lgth5, ncol=lgth2)

for (i in 1:lgth2){
  div[,i]<-(sp_cols[,i]/dat1$effort)*500
  div[,i]<-round(div[,i],0)
}
div<-as.data.frame(div)
colnames(div)<-caps$sp

for(i in 1:lgth2){
  div[,i]<-as.integer(div[,i])
}
dat4<-data.frame(dat1[,1:8],div,dat1[(length(dat1)-3):length(dat1)])
head(dat1, 3)
head(dat4, 3)

compare<-data.frame(site=dat1$obs,raw=dat1$Pog,eff=dat1$effort,std=dat4$Pog)
head(compare)

# plot(compare$raw, compare$std)
# cor.test(compare$raw, compare$std)


# THE DATA:

# The mod_type can also be used to specify which data set to use:

# 1 = Hdat (Hincks only, all sites)
# 2 = Pdat (Pinks only, all sites)
# 3 = dat1 (both locations, all sites)

# 1:
Hdat<-dat4[dat4$location=="Hincks",]
Hdat<-droplevels(Hdat)
rownames(Hdat)<-1:nrow(Hdat)
Hdat$obs_effect<-as.numeric(1:nrow(Hdat))
head(Hdat)

# 2:
Pdat<-dat4[dat4$location=="Pinkawillinie",]
Pdat<-droplevels(Pdat)
rownames(Pdat)<-1:nrow(Pdat)
Pdat$obs_effect<-as.numeric(1:nrow(Pdat))
head(Pdat)

# 3: dat4
dat4$obs_effect<-as.numeric(1:nrow(dat4))
head(dat4)


head(Hdat)
levels(Hdat$site)

newH<-Hdat[-which(Hdat$site %in% c("I1", "I2", "I3", "I4", "I7")),]
newH<-droplevels(newH)
levels(newH$site)
head(newH); dim(newH)

comb_dat<-as.data.frame(apply(newH[,c(which(colnames(newH)=="A_ina"):which(colnames(newH)=="D_dam"))],2,function(x) aggregate(x~newH[,which(colnames(newH)=="site")],FUN=sum)))
comb_dat<-comb_dat[,grep(".x", colnames(comb_dat))]
colnames(comb_dat)<-colnames(newH)[which(colnames(newH)=="A_ina"):which(colnames(newH)=="D_dam")]
comb_dat<-comb_dat[,1:which(colnames(comb_dat)=="Pog")]
comb_dat<-data.frame(site=unique(as.character(newH[,which(colnames(newH)=="site")])),comb_dat)

# total abundance:
comb_dat$abund<-rowSums(comb_dat[,2:ncol(comb_dat)])

# species richness:
comb_dat$sp_rich<-rowSums(t(apply(comb_dat[,c(which(colnames(comb_dat)=="A_ina"):which(colnames(comb_dat)=="Pog"))],1,function(x) ifelse(x==0,0,1)))) 
head(comb_dat); dim(comb_dat)

# relative abundance:
rel_abund<-data.frame(site=comb_dat$site,round(apply(comb_dat[,c(which(colnames(comb_dat)=="A_ina"):which(colnames(comb_dat)=="Pog"))],2,function(x) x/comb_dat$abund),3), abund=comb_dat$abund, sp_rich=comb_dat$sp_rich)

# simpson's index:
rel_abund$simps_ind<-1-(rowSums(rel_abund[,c(which(colnames(rel_abund)=="A_ina"):which(colnames(rel_abund)=="Pog"))]^2))

# reciprocal simpson's index:
rel_abund$simps_recip<-1/rel_abund$simps_ind
head(rel_abund); dim(rel_abund)

library(vegan)
diversity(rel_abund[,c(which(colnames(rel_abund)=="A_ina"):which(colnames(rel_abund)=="Pog"))], index="simpson")

# barplot(rel_abund$simps_ind[c(3:6,1:2)])


# full data set:

dat2<-read.table("all_reptile_data.txt",header=T)
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

eft<-data.frame(aggregate(xdat$effort~xdat$site, FUN=sum))
colnames(eft)<-c("site","effort")
eft

head(sp_div,3); dim(sp_div)
rownames(sp_div)==as.character(eft$site)

sdiv<-sp_div

# captures / 1000 trap nights
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


### simpson's index:

quartz("", 8,8, dpi=80, pointsize=16)
par(mfrow=c(2,2), mar=c(4,4,0.5,0.5), mgp=c(2.9,1,0))
plot(sum_dat$fire_cat, sum_dat$sp_rich, xlab="", las=1, ylab="Number of species")
plot(sum_dat$fire_cat, sum_dat$simps_ind2, xlab="", las=1, ylab="Simpson's Index")

plot(sum_dat$fire_cat, sum_dat$even, xlab="", las=1, ylab="Evenness")
plot(sum_dat$simps_ind2, sum_dat$even, xlab="", las=1, ylab="Evenness", pch=20)
title(xlab="Simpson's Index",mgp=c(2.5,1,0))

### shannon's index:

quartz("", 8,8, dpi=80, pointsize=16)
par(mfrow=c(2,2), mar=c(4,4,0.5,0.5), mgp=c(2.9,1,0))

plot(sum_dat$fire_cat, sum_dat$sp_rich, xlab="", las=1, ylab="Number of species")

plot(sum_dat$fire_cat, sum_dat$shann_ind, xlab="", las=1, ylab="Shannon's Index")

plot(sum_dat$fire_cat, sum_dat$even2, xlab="", las=1, ylab="Evenness")

plot(sum_dat$shann_ind, sum_dat$even2, xlab="", las=1, ylab="Evenness", pch=20)
title(xlab="Shannon's Index",mgp=c(2.5,1,0))

# compare simpson and shannon:

quartz("", 6,6, dpi=70, pointsize=20)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(sum_dat$simps_ind2, sum_dat$shann_ind, xlab="", las=1, ylab="Shannon's Index", pch=20)
text(4, 2.4,paste("r = ",round(cor.test(sum_dat$simps_ind2, sum_dat$shann_ind)$estimate,2),sep=""), adj=0)
title(xlab="Simpson's Index",mgp=c(2.5,1,0))


# simpson's index should be less related to richness than shannon:
# shannon index puts more weight on species richness - rare species count more
# simpson's index puts more weight on dominance - dominant species count more
cor.test(sum_dat$simps_ind2, sum_dat$sp_rich)
cor.test(sum_dat$shann_ind, sum_dat$sp_rich)

quartz("", 8,4, dpi=80, pointsize=16)
par(mfrow=c(1,2), mar=c(4,4,0.5,0.5), mgp=c(2.9,1,0))

plot(sum_dat$sp_rich, sum_dat$simps_ind2)
plot(sum_dat$sp_rich, sum_dat$shann_ind)



# but since the simpson index is a measure of dominance, should it have a stronger relationship (confounding?) with evenness? No, they're the same:
cor.test(sum_dat$simps_ind2, sum_dat$even)
cor.test(sum_dat$shann_ind, sum_dat$even2)



# Rank abundance:
# plot(sum_dat$simps_ind, sum_dat$shann_ind)

head(rel_abund,3); dim(rel_abund)
head(sp_div2, 3); dim(sp_div2)

sp_abund<-data.frame(species=names(colSums(sp_div2)), total_abund=colSums(sp_div2))
sp_abund<-sp_abund[order(sp_abund$total_abund),]
rownames(sp_abund)<-1:nrow(sp_abund)
head(sp_abund); dim(sp_abund)

quartz("",12,4,pointsize=16)
par(mar=c(4,4,1,1), mgp=c(2.5,1,0))
hist(sp_abund$total_abund, breaks=seq(0,500,by=10), main="", ylab="Number of species", xlab="",las=1, col="grey80")
title(xlab="Abundance", mgp=c(2.5,1,0))


# how many rare species in each fire category?
rare_sp<-as.character(sp_abund$species[sp_abund$total_abund<5])
sp_div3<-sp_div2[,rare_sp]

sum_dat$rare_sp<-apply(sp_div3, 1, function(x) length(which(x>0)))

quartz("", 6,6, dpi=70, pointsize=20)
par(mfrow=c(1,1), mar=c(3,4,1,1), mgp=c(2.5,1,0))
plot(sum_dat$fire_cat, sum_dat$rare_sp, xlab="", ylab="Number of rare species (1-4 captures)", las=1)

# plot(sum_dat$simps_ind2, sum_dat$rare_sp)
# plot(sum_dat$shann_ind, sum_dat$rare_sp)

cor.test(sum_dat$shann_ind, sum_dat$rare_sp)
cor.test(sum_dat$simps_ind2, sum_dat$rare_sp)

colnames(sp_div2)


# how many dominant (high abundance) species in each fire category?
dom_sp<-as.character(sp_abund$species[sp_abund$total_abund>200])
sp_div4<-sp_div2[,dom_sp]

sp_div5<-data.frame(sp_div4)
sp_div5$site<-rownames(sp_div4)
rownames(sp_div5)<-1:nrow(sp_div5)

sdat<-data.frame(site=sum_dat$site,fire_cat=sum_dat$fire_cat)

sp_div5<-merge(sp_div5, sdat)


quartz("", 8,8, dpi=80, pointsize=16)
par(mfrow=c(2,2), mar=c(4,4,0.5,0.5), mgp=c(2.9,1,0))
plot(sp_div5$fire_cat, sp_div5$C_euc, las=1, xlab="", ylab=expression(italic("Ctenotus euclae")))
plot(sp_div5$fire_cat, sp_div5$E_in, las=1, xlab="", ylab=expression(italic("Egernia inornata")))
plot(sp_div5$fire_cat, sp_div5$N_stel, las=1, xlab="", ylab=expression(italic("Nephrurus stellatus")))
plot(sp_div5$fire_cat, sp_div5$N_stel, las=1, xlab="", ylab=expression(italic("Ctenophorus fordi")))


# analyse

head(sum_dat)
sum_dat$location<-NA
sum_dat$location[grep("H",sum_dat$site)]<-"Hincks"
sum_dat$location[grep("P",sum_dat$site)]<-"Pinks"

m1<-lm(shann_ind~fire_cat+location, data=sum_dat)
summary(m1)
anova(m1)

library(lme4)

c1<-matrix(data=c(rep(5,4),c(16,1,1,2)), nrow=2, byrow=T)
c1
diversity(c1, index="shannon")





