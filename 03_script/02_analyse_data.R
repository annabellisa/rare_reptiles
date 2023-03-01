

# 02_analyse_data

load("04_workspaces/processed_data.RData")

# processed data of 14 sites with species diversity metrics
head(sum_dat, 6);dim(sum_dat)

### simpson's index:

dev.new(height = 8, width = 8, noRStudioGD = T, dpi=80, pointsize=16)
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





