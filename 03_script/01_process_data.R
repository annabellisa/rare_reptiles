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

library("lme4"); library("vegan")

# Set up the data:

dat1<-read.table("01_data/sp_summaries.txt",header=T)
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

diversity(rel_abund[,c(which(colnames(rel_abund)=="A_ina"):which(colnames(rel_abund)=="Pog"))], index="simpson")

# barplot(rel_abund$simps_ind[c(3:6,1:2)])

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

dir()
save.image("04_workspaces/processed_data.RData")


