# ------------------------------------ #
# ---------- RARE REPTILES  ---------- #
# ------------------------------------ #

### Analysis fire effects on rarity and dominance in a reptile community 

# ---------- 02_rarefaction  ---------- #

### Script authors: Amber Lim & Annabel Smith

# Load and tidy workspace and remove everything except necessary objects:
# load("04_workspaces/analysed_data.RData"); rm(list=setdiff(ls(), c("sum_dat","rare.data","sp_div2","sp_div_25","sp_div_5")))

# load workspace
load("04_workspaces/rarefaction.RData")

# Load functions:
invisible(lapply(paste("02_functions/",dir("02_functions"),sep=""), function(x) source(x)))

# Following Hsei et al. 2022, iNEXT introduction
library(devtools)
# install_github('AnneChao/iNEXT')

## import packages
library(iNEXT)
library(ggplot2)

# save.image("04_workspaces/rarefaction.RData")

### DATA ### 

# ----

# full processed diversity data:
head(sum_dat, 6);dim(sum_dat)

# rare species richness & abundance
head(rare.data); dim(rare.data)

# Full site x species matrix (14 x 41)
head(sp_div2, 3); dim(sp_div2)

# site x species matrix, rare 25% species
head(sp_div_25, 3); dim(sp_div_25)

# site x species matrix, rare 5% species
head(sp_div_5, 3); dim(sp_div_5)

# ----

# Analysis

# Full site x species matrix (14 x 41)
head(sum_dat,3)
levels(sum_dat$fire_cat)

# Aggregate abundances across burn categories (equivalent to the bird data in Hsing et al. 2022):
b_sites<-as.character(sum_dat$site[sum_dat$fire_cat=="Burnt"])
m_sites<-as.character(sum_dat$site[sum_dat$fire_cat=="Medium"])
ub_sites<-as.character(sum_dat$site[sum_dat$fire_cat=="Unburnt"])

head(sp_div2, 3); dim(sp_div2)

abund_mat<-data.frame(Burnt=rowSums(t(sp_div2[which(rownames(sp_div2) %in% b_sites),])))
abund_mat$Medium<-rowSums(t(sp_div2[which(rownames(sp_div2) %in% m_sites),]))
abund_mat$Unburnt<-rowSums(t(sp_div2[which(rownames(sp_div2) %in% ub_sites),]))

# aggregated abundance data:
head(abund_mat); dim(abund_mat)

# run iNEXT analysis:
reptile_iN<-iNEXT(abund_mat, q=0, datatype="abundance")
str(reptile_iN)
r_dat<-reptile_iN$DataInfo
est_dat<-reptile_iN$iNextEst
asy_dat<-reptile_iN$AsyEst
str(est_dat)

head(asy_dat)

asy_dat$Assemblage<-factor(asy_dat$Assemblage,levels=c("Unburnt","Medium","Burnt"))

dev.new(width=8, height=5, dpi=80, pointsize=20, noRStudioGD = T)
par(mfrow=c(1,1), mar=c(2.5,4,0,1), mgp=c(2.8,0.8,0), oma=c(0,0,1,7))
  
  # species richness (5% of max.)
  plot(c(1:3),asy_dat$Estimator[asy_dat$Diversity=="Shannon diversity"], xlim=c(0.5,3.5), pch=20, col="grey50", xaxt="n",ylim= c((min(asy_dat$LCL)),max(asy_dat$UCL)),ylab="Asymptotic diversity estimate",xlab="", las = 1, cex = 1)
  
  arrows(c(1:3),asy_dat$LCL[asy_dat$Diversity=="Shannon diversity"],c(1:3),asy_dat$UCL[asy_dat$Diversity=="Shannon diversity"],length=0.03,code=3,angle=90)
  
  points(c(1:3)-0.2,asy_dat$Estimator[asy_dat$Diversity=="Species richness"], pch=20, col="black")
  arrows(c(1:3)-0.2,asy_dat$LCL[asy_dat$Diversity=="Species richness"],c(1:3)-0.2,asy_dat$UCL[asy_dat$Diversity=="Species richness"],length=0.03,code=3,angle=90)
  
  points(c(1:3)+0.2,asy_dat$Estimator[asy_dat$Diversity=="Simpson diversity"], pch=20, col="grey75")
  arrows(c(1:3)+0.2,asy_dat$LCL[asy_dat$Diversity=="Simpson diversity"],c(1:3)+0.2,asy_dat$UCL[asy_dat$Diversity=="Simpson diversity"],length=0.03,code=3,angle=90)
  
  axis(1,at=c(1:3),labels=F)
  axis(1,at=c(1:3), cex.axis=1,labels=levels(asy_dat$Assemblage),tick=F)
  
  # plot legend:
par(xpd=NA)
legend(x=3.75,y=max(asy_dat$UCL)+1.1, title = "Hill number order", legend = c(expression(paste(italic("q")," = 0", sep="")), expression(paste(italic("q")," = 1",sep="")),"q = 2"), pt.cex = 1, pch = c(20), bty = "n",adj=0, title.adj=0, col=c("black","grey50","grey75"))
par(xpd=F)








# iNEXT tutorial examples:

# abundance data for two treatments, no species identity
data(spider)
str(spider)
head(spider)
spider_iN<-iNEXT(spider, q=0, datatype="abundance")
spider_iN$DataInfo

# abundance data for 41 species at two locations:
data(bird)
str(bird)
head(bird); dim(bird)
bird_iN<-iNEXT(bird, q=0, datatype="abundance")

# frequency distribution for five sampling locations, no species identity
data(ant)
str(ant)
head(ant)
ant_iN<-iNEXT(ant, q=0, datatype="incidence_freq")

# presence absence data (0,1) for three locations, species x sampling location
data(ciliates)
str(ciliates)
head(ciliates)
head(ciliates$EtoshaPan); dim(ciliates$EtoshaPan)
head(ciliates$CentralNamibDesert); dim(ciliates$CentralNamibDesert)
range(ciliates$EtoshaPan)
range(ciliates$CentralNamibDesert)
ciliate_iN<-iNEXT(ciliates, q=0, datatype="incidence_raw")


