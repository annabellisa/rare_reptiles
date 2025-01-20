# ------------------------------------ #
# ---------- RARE REPTILES  ---------- #
# ------------------------------------ #

### Analysis fire effects on rarity and dominance in a reptile community 

# ---------- 03_rarefaction  ---------- #

### Script authors: Annabel Smith & Amber Lim

# Load and tidy workspace and remove everything except necessary objects:
# load("04_workspaces/analysed_data.RData"); rm(list=setdiff(ls(), c("sum_dat","rare.data","sp_div2","sp_div_25","sp_div_5")))

# load workspace
load("04_workspaces/rarefaction.RData")

# Load functions:
invisible(lapply(paste("02_functions/",dir("02_functions"),sep=""), function(x) source(x)))

library("RColorBrewer")

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

# Analysis ----

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
reptile_iN<-iNEXT(abund_mat, q=c(0,1,2), datatype="abundance",nboot=200, knots=120)

# save.image("04_workspaces/rarefaction.RData")

str(reptile_iN)

r_dat<-reptile_iN$DataInfo
est_dat<-reptile_iN$iNextEst
asy_dat<-reptile_iN$AsyEst

# DataInfo: 
# n = reference sample size, S.obs = observed species richness, SC = sample coverage estimate for the reference sample, fi = first ten frequency counts
r_dat

# iNextEst:
# Diversity estimates and related statistics

# there are size based and coverage based estimates, both containing CIs obtained from a bootstrap method. Sampling uncertainty is greater in the coverage based estimates, so CIs are wider for this group of estimates. 

# m = sample size based on 40 equally spaced knots, from 1 to (2 x n), where n is the reference sample size for each group. This locates the reference sample at the mid point of the diversity estimates. qD is the Hill number diversity estimate of order q; SC is the sample coverage. Both of these estimates have associated CIs. 

seq(1,(558*2), length.out=40)
seq(1,(168*2), length.out=40)

str(est_dat)

head(est_dat$size_based,3); dim(est_dat$size_based)
head(est_dat$coverage_based,3); dim(est_dat$coverage_based)

dim(est_dat$size_based[which(est_dat$size_based=="Burnt"),])

# AsyEst
# Asymptotic diversity estimates for richness, Shannon and Simpson diversity
head(asy_dat,3); dim(asy_dat)

# save.image("04_workspaces/rarefaction.RData")

# ----

# Plot iNext object ----

# Plot size-based, then coverage-based
e_size<-est_dat$size_based
e_cov<-est_dat$coverage_based
head(e_cov)
head(e_size,3); dim(e_size)

e_size$Assemblage<-factor(e_size$Assemblage,levels=c("Unburnt","Medium","Burnt"))
e_cov$Assemblage<-factor(e_cov$Assemblage,levels=c("Unburnt","Medium","Burnt"))

head(e_size,3); dim(e_size)
head(e_cov,3); dim(e_cov)

e_size[which(e_size$Method=="Observed"),]
e_cov[which(e_cov$Method=="Observed"),]

e_size[e_size$m<100,]

# q.col.test<-hcl.colors(16, palette = "viridis", alpha = 1)
# plot(1:16, 1:16, pch=20, col=q.col.test, cex=3)

# Plot by order, with each assemblage shown on the same panel, with different colours
orders.q<-unique(e_size$Order.q)
assemb.q<-unique(e_size$Assemblage)

# Set colours:
q.col<-hcl.colors(16, palette = "viridis", alpha = 0.5)[c(1,7,14)]
assemb.col<-hcl.colors(16, palette = "viridis", alpha = 1)[c(1,7,14)]

# Plot
dev.new(width=10.5, height=7, dpi=80, pointsize=18, noRStudioGD = T)
par(mfrow=c(2,3), mar=c(4.5,4,1,1), mgp=c(2,0.8,0), oma=c(0,0,1,7))

# sample-size-based

for(i in 1:length(orders.q)){

order.thisrun<-orders.q[i]
data.thisrun<-e_size[e_size$Order.q==order.thisrun,]

# Chao et al. 2014 allow the ylim to vary for each q, whereas Hsieh et al. 2016 keep them consistent. The first of these sets a consistent ylim, as the q with the maximum number of species. The second allows it to vary for each q. The third sets the max ylim at 20 to enable interpretation. 
y.lim.all<-c(min(e_size$qD.LCL),max(e_size$qD.UCL))
y.lim.ind<-c(min(data.thisrun$qD.LCL),max(data.thisrun$qD.UCL))
if(i==1) y.lim.var<-c(min(e_size$qD.LCL),max(e_size$qD.UCL)) else y.lim.var<-c(min(e_size$qD.LCL),20)

plot(data.thisrun$m[data.thisrun$Assemblage=="Burnt"],data.thisrun$qD[data.thisrun$Assemblage=="Burnt"], type="n", ylim=y.lim.var,xlim=c(0,1500), las=1, xlab="Number of individuals",ylab="Estimated species diversity")

pg.ci(x="m","data.thisrun",x.subset="Assemblage",colour=q.col, lower="qD.LCL",upper="qD.UCL")

mtext(bquote("("*.(letters[i])*") Sample-size-based, "*italic("q")*" = "*.(order.thisrun)), adj=0, cex=0.625)
                        
for (j in 1:length(assemb.q)){

  asb.thisrun<-levels(assemb.q)[j]
  assemb.thisrun<-data.thisrun[data.thisrun$Assemblage==asb.thisrun,]
  lines(assemb.thisrun$m[assemb.thisrun$Method=="Rarefaction"],assemb.thisrun$qD[assemb.thisrun$Method=="Rarefaction"])
  lines(assemb.thisrun$m[assemb.thisrun$Method=="Extrapolation"],assemb.thisrun$qD[assemb.thisrun$Method=="Extrapolation"], lty=2)
  points(assemb.thisrun$m[assemb.thisrun$Method=="Observed"],assemb.thisrun$qD[assemb.thisrun$Method=="Observed"], pch=20,col=assemb.col[j], cex=2)
  
  # add.max<- ifelse(i==1, 40, 20) 
  # axis(side=2,at=add.max,labels=add.max,tick=F,font=2, las=1)
  
} # close j assemblage

} # close order size-based

# plot legend:
par(xpd=NA)
legend(x=1600,y=max(y.lim.var), title = "Fire category", legend = c("Long","Medium","Recent"), pt.cex = 2, pch = c(20, 20, 20), bty = "n", title.adj=0,col=assemb.col)
legend(x=1550,y=max(y.lim.var)-7, title = "", legend = unique(e_size$Method)[c(1,3)], lty = c(1,2), bty = "n", title.adj=0)
par(xpd=F)

# coverage-based

for(i in 1:length(orders.q)){
  
  order.thisrun<-orders.q[i]
  data.thisrun<-e_cov[e_cov$Order.q==order.thisrun,]
  
  # Chao et al. 2014 allow the ylim to vary for each q, whereas Hsieh et al. 2016 keep them consistent. The first of these sets a consistent ylim, as the q with the maximum number of species. The second allows it to vary for each q. The third sets the max ylim at 20 to enable interpretation. 
  y.lim.all<-c(min(e_size$qD.LCL),max(e_size$qD.UCL))
  y.lim.ind<-c(min(data.thisrun$qD.LCL),max(data.thisrun$qD.UCL))
  if(i==1) y.lim.var<-c(min(e_size$qD.LCL),max(e_size$qD.UCL)) else y.lim.var<-c(min(e_size$qD.LCL),20)
  
  plot(data.thisrun$SC[data.thisrun$Assemblage=="Burnt"],data.thisrun$qD[data.thisrun$Assemblage=="Burnt"], type="n", ylim=y.lim.var,xlim=c(min(e_cov$SC),max(e_cov$SC)), las=1, xlab="Sample coverage",ylab="Estimated species diversity")
  
  pg.ci(x="SC","data.thisrun",x.subset="Assemblage",colour=q.col, lower="qD.LCL",upper="qD.UCL")
  
  mtext(bquote("("*.(letters[i+3])*") Coverage-based, "*italic("q")*" = "*.(order.thisrun)), adj=0, cex=0.625)
  
  for (j in 1:length(assemb.q)){
    
    asb.thisrun<-levels(assemb.q)[j]
    assemb.thisrun<-data.thisrun[data.thisrun$Assemblage==asb.thisrun,]
    lines(assemb.thisrun$SC[assemb.thisrun$Method=="Rarefaction"],assemb.thisrun$qD[assemb.thisrun$Method=="Rarefaction"])
    lines(assemb.thisrun$SC[assemb.thisrun$Method=="Extrapolation"],assemb.thisrun$qD[assemb.thisrun$Method=="Extrapolation"], lty=2)
    points(assemb.thisrun$SC[assemb.thisrun$Method=="Observed"],assemb.thisrun$qD[assemb.thisrun$Method=="Observed"], pch=20,col=assemb.col[j], cex=2)
    
  } # close j assemblage
  
} # close order

# ----

# Plot sample coverage against the number of individuals 

# From Chao et al. 2014: The sample-size- and coverage-based curves are linked by a sample completeness curve (Figs. 4 and 7), which reveals the relationship between sample size (number of individuals or number of sampling units) and sample completeness. This curve illustrates how much sampling effort is needed to achieve a pre-determined level of sample completeness.

# save.image("04_workspaces/rarefaction.RData")

dev.new(width=12, height=5, dpi=80, pointsize=18, noRStudioGD = T)
par(mfrow=c(1,2), mar=c(4.5,4,1,1), mgp=c(2.3,0.8,0), oma=c(0,0,1,7))

# All q orders have the same sample coverage, so we only need a single plot to show the sample coverage

i=1
  
  order.thisrun<-orders.q[i]
  data.thisrun<-e_size[e_size$Order.q==order.thisrun,]
  head(data.thisrun,2)
  
  plot(data.thisrun$m[data.thisrun$Assemblage=="Burnt"],data.thisrun$SC[data.thisrun$Assemblage=="Burnt"], type="n", ylim=c(min(data.thisrun$SC.LCL), max(data.thisrun$SC.UCL)),xlim=c(1,1500), las=1, xlab="Number of individuals",ylab="Sample coverage")
  
  pg.ci(x="m","data.thisrun",x.subset="Assemblage",colour=q.col,lower="SC.LCL",upper="SC.UCL")
  
  for (j in 1:length(assemb.q)){
    
    asb.thisrun<-levels(assemb.q)[j]
    assemb.thisrun<-data.thisrun[data.thisrun$Assemblage==asb.thisrun,]
    lines(assemb.thisrun$m[assemb.thisrun$Method=="Rarefaction"],assemb.thisrun$SC[assemb.thisrun$Method=="Rarefaction"])
    lines(assemb.thisrun$m[assemb.thisrun$Method=="Extrapolation"],assemb.thisrun$SC[assemb.thisrun$Method=="Extrapolation"], lty=2)
    
  } # close j assemblage
  
  # plot points separately so lines don't overlap:
  
  j=1
  
  for (j in 1:length(assemb.q)){
  
    asb.thisrun<-levels(assemb.q)[j]
  assemb.thisrun<-data.thisrun[data.thisrun$Assemblage==asb.thisrun,]
  points(assemb.thisrun$m[assemb.thisrun$Method=="Observed"],assemb.thisrun$SC[assemb.thisrun$Method=="Observed"], pch=20,col=assemb.col[j], cex=2)
    
  } # close j assemblage
  
  mtext(paste("(",(letters[i]),")",sep=""), adj=0, cex=1, line=0.4)
  
  # Plot another one, zoomed in on the x-axis:
  
  i=1
  
  order.thisrun<-orders.q[i]
  data.thisrun<-e_size[e_size$Order.q==order.thisrun,]
  head(data.thisrun,2)
  
  plot(data.thisrun$m[data.thisrun$Assemblage=="Burnt"],data.thisrun$SC[data.thisrun$Assemblage=="Burnt"], type="n", ylim=c(min(data.thisrun$SC.LCL), max(data.thisrun$SC.UCL)),xlim=c(1,300), las=1, xlab="Number of individuals",ylab="Sample coverage", xaxt="n")
  axis(side=1, at=c(0,100,200,300), labels=T)
  
  pg.ci(x="m","data.thisrun",x.subset="Assemblage",colour=q.col,lower="SC.LCL",upper="SC.UCL")
  
  for (j in 1:length(assemb.q)){
    
    asb.thisrun<-levels(assemb.q)[j]
    assemb.thisrun<-data.thisrun[data.thisrun$Assemblage==asb.thisrun,]
    lines(assemb.thisrun$m[assemb.thisrun$Method=="Rarefaction"],assemb.thisrun$SC[assemb.thisrun$Method=="Rarefaction"])
    lines(assemb.thisrun$m[assemb.thisrun$Method=="Extrapolation"],assemb.thisrun$SC[assemb.thisrun$Method=="Extrapolation"], lty=2)
    points(assemb.thisrun$m[assemb.thisrun$Method=="Observed"],assemb.thisrun$SC[assemb.thisrun$Method=="Observed"], pch=20,col=assemb.col[j], cex=2)
    
  } # close j assemblage
  
  mtext(paste("(",(letters[i+1]),")",sep=""), adj=0, cex=1, line=0.4)
  
# plot legend:
par(xpd=NA)
legend(x=320,y=1, title = "Fire category", legend = c("Long","Medium","Recent"), pt.cex = 2, pch = c(20, 20, 20), bty = "n", title.adj=0,col=assemb.col)
legend(x=300,y=0.6, title = "", legend = unique(e_size$Method)[c(1,3)], lty = c(1,2), bty = "n", title.adj=0)
par(xpd=F)

# save.image("04_workspaces/rarefaction.RData")





# Run ggINEXT to check custom plots against inbuilt plotting functions:
r_iN1<-iNEXT(abund_mat, q=c(0,1,2), datatype="abundance")
est_iN1<-r_iN1$iNextEst
head(est_iN1$size_based,3); dim(est_iN1$size_based)

est_iN1$size_based[which(est_iN1$size_based=="Burnt"),1:4]

ggiNEXT(r_iN1, type=1, se=TRUE, facet.var="Assemblage", color.var="Order.q", grey=FALSE)
ggiNEXT(r_iN1, type=1, se=TRUE, facet.var="Order.q", color.var="Assemblage", grey=FALSE)
ggiNEXT(r_iN1, type=3, se=TRUE, facet.var="Order.q", color.var="Assemblage", grey=FALSE)
ggiNEXT(r_iN1, type=3, se=TRUE, facet.var="Assemblage", color.var="Order.q", grey=FALSE)
ggiNEXT(r_iN1, type=2, se=TRUE, facet.var="None", color.var="Assemblage", grey=FALSE)

# Plot asymptotic diversity estimate:

asy_dat$Assemblage<-factor(asy_dat$Assemblage,levels=c("Unburnt","Medium","Burnt"))

dev.new(width=6, height=4, dpi=80, pointsize=16, noRStudioGD = T)
par(mfrow=c(1,1), mar=c(2.5,4,0,1), mgp=c(2.8,0.8,0), oma=c(0,0,1,7))
  
  # species richness (5% of max.)
  plot(c(1:3),asy_dat$Estimator[asy_dat$Diversity=="Shannon diversity"], xlim=c(0.5,3.5), pch=20, col="grey50", xaxt="n",ylim= c((min(asy_dat$LCL)),max(asy_dat$UCL)),ylab="Asymptotic diversity estimate",xlab="", las = 1, cex = 1)
  
  arrows(c(1:3),asy_dat$LCL[asy_dat$Diversity=="Shannon diversity"],c(1:3),asy_dat$UCL[asy_dat$Diversity=="Shannon diversity"],length=0.03,code=3,angle=90)
  
  points(c(1:3)-0.2,asy_dat$Estimator[asy_dat$Diversity=="Species richness"], pch=20, col="black")
  arrows(c(1:3)-0.2,asy_dat$LCL[asy_dat$Diversity=="Species richness"],c(1:3)-0.2,asy_dat$UCL[asy_dat$Diversity=="Species richness"],length=0.03,code=3,angle=90)
  
  points(c(1:3)+0.2,asy_dat$Estimator[asy_dat$Diversity=="Simpson diversity"], pch=20, col="grey75")
  arrows(c(1:3)+0.2,asy_dat$LCL[asy_dat$Diversity=="Simpson diversity"],c(1:3)+0.2,asy_dat$UCL[asy_dat$Diversity=="Simpson diversity"],length=0.03,code=3,angle=90)
  
  axis(1,at=c(1:3),labels=F)
  axis(1,at=c(0.8,2,3.2), cex.axis=1,labels=levels(asy_dat$Assemblage),tick=F)
  
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


