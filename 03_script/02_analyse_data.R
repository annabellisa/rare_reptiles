

# 02_analyse_data

load("04_workspaces/processed_data.RData")

library("MASS")

# species matrix into data frame for processing
sp_div2 <- as.data.frame(sp_div2)
#save.image("04_workspaces/processed_data.RData")

# captures / 1000 trap nights. 2 seasons of data combined
head(sp_div2, 3); dim(sp_div2)

# rearranging so unburnt is baseline
sum_dat$fire_cat <- factor(sum_dat$fire_cat, levels = c("Unburnt", "Medium", "Burnt"))

# making site a factor
sum_dat$site <- as.factor(sum_dat$site)

#List of 14 sites, taking site column and calling "H" to create 2 different variables for site, either "H" or "P"
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

sp_5 <- sp_sum2$Species_Code[which(sp_sum2$max_5 == "r")]

#pulling out data set for rare 25% species
sp_div_25 <- sp_div2[,which(colnames(sp_div2)%in%sp_25)]
head(sp_div_25, 3); dim(sp_div_25)

#pulling out data set for rare 5% species
sp_div_5 <- sp_div2[,which(colnames(sp_div2)%in%sp_5)]
head(sp_div_5, 3); dim(sp_div_5)

rare.data <- data.frame(site=names(rowSums(sp_div_25)), abund_25 = rowSums(sp_div_25), pres_25 = ifelse(rowSums(sp_div_25) == 0, 0, 1), abund_5 = rowSums(sp_div_5), pres_5 = ifelse(rowSums(sp_div_5) == 0, 0, 1))

row.names(rare.data) <- 1:nrow(rare.data)
head(sum_dat, 6);dim(sum_dat)

#combining data sets into sum_dat
sum_dat <- merge(sum_dat, rare.data, by="site", all.x = T, all.y = F)

#data fully processed and ready to process 9th May 2023 in sum_dat

# species richness as a function of fire category (generalized linear model with negative binomial linear structure). Statistical model = no negative results because count data. 
#Effect of fire can vary depending on location - fire cat*location have no effect on species richness, fire cat AND location have no effect on species richness.
m1_a <- glm.nb(sp_rich~fire_cat*location,data = sum_dat)
m1_b <- glm.nb(sp_rich~fire_cat+location,data = sum_dat)
m1_c <- glm.nb(sp_rich~fire_cat,data = sum_dat)

summary(m1_a);anova(m1_a)
summary(m1_b);anova(m1_b)
summary(m1_c);anova(m1_c)

#save.image("04_workspaces/processed_data.RData")

head(sum_dat, 6);dim(sum_dat)

summary(sum_dat$simps_ind2) 

#do glm.nb for abund_25 and abund_5
#do glm for shann_ind
#do glm for pres_25 and pres_5, change "Gamma" to "binomial"

#species diversity as simpson's index:
m2_a <- glm(simps_ind2~fire_cat*location,data = sum_dat, family = "Gamma")
m2_b <- glm(simps_ind2~fire_cat+location,data = sum_dat, family = "Gamma")
m2_c <- glm(simps_ind2~fire_cat,data = sum_dat, family = "Gamma")
m2_d <- glm(simps_ind2~1,data = sum_dat, family = "Gamma")

#Likelihood ratio test - compares more complicated model (interaction model) with simpler model 
anova(m2_a, m2_b, test = "F") # p value for interaction term = 0.55 (fire x location)
anova(m2_b, m2_c, test = "F") # p value for location term = 0.10
anova(m2_c, m2_d, test = "F") # p value for fire term = 0.06

summary(m2_a);anova(m2_a)
summary(m2_b);anova(m2_b)
summary(m2_c);anova(m2_c)
summary(m2_d);anova(m2_d)

summary(sum_dat$shann_ind)

### simpson's index:

dev.new(height = 8, width = 8, noRStudioGD = T, dpi=80, pointsize=16)
par(mfrow=c(2,2), mar=c(4,4,0.5,0.5), mgp=c(2.9,1,0))
plot(sum_dat$fire_cat, sum_dat$sp_rich, xlab="", las=1, ylab="Number of species")
plot(sum_dat$fire_cat, sum_dat$simps_ind2, xlab="", las=1, ylab="Simpson's Index")

plot(sum_dat$fire_cat, sum_dat$even, xlab="", las=1, ylab="Evenness")
plot(sum_dat$simps_ind2, sum_dat$even, xlab="", las=1, ylab="Evenness", pch=20)
title(xlab="Simpson's Index",mgp=c(2.5,1,0))

### shannon's index:
#quartz("", 8,8, dpi=80, pointsize=16)
dev.new(height = 8, width = 8, noRStudioGD = T, dpi=80, pointsize=16)
par(mfrow=c(2,2), mar=c(4,4,0.5,0.5), mgp=c(2.9,1,0))

plot(sum_dat$fire_cat, sum_dat$sp_rich, xlab="", las=1, ylab="Number of species")

plot(sum_dat$fire_cat, sum_dat$shann_ind, xlab="", las=1, ylab="Shannon's Index")

plot(sum_dat$fire_cat, sum_dat$even2, xlab="", las=1, ylab="Evenness")

plot(sum_dat$shann_ind, sum_dat$even2, xlab="", las=1, ylab="Evenness", pch=20)
title(xlab="Shannon's Index",mgp=c(2.5,1,0))

head(sum_dat, 6);dim(sum_dat)
str(sum_dat)

# determining type of data for different statistical analysis types
# sp_rich data = count; appropriate model is negative binomial distribution
summary(sum_dat$sp_rich)

summary(sum_dat$simps_ind2)

summary(sum_dat$even)

summary(sum_dat$shann_ind)

# compare simpson and shannon:
#quartz("", 6,6, dpi=70, pointsize=20)

dev.new(height = 8, width = 8, noRStudioGD = T, dpi=80, pointsize=16)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(sum_dat$simps_ind2, sum_dat$shann_ind, xlab="", las=1, ylab="Shannon's Index", pch=20)
text(4, 2.4,paste("r = ",round(cor.test(sum_dat$simps_ind2, sum_dat$shann_ind)$estimate,2),sep=""), adj=0)
title(xlab="Simpson's Index",mgp=c(2.5,1,0))


# simpson's index should be less related to richness than shannon:
# shannon index puts more weight on species richness - rare species count more
# simpson's index puts more weight on dominance - dominant species count more
cor.test(sum_dat$simps_ind2, sum_dat$sp_rich)
cor.test(sum_dat$shann_ind, sum_dat$sp_rich)

dev.new(height = 8, width = 8, noRStudioGD = T, dpi=80, pointsize=16)
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

#quartz("",12,4,pointsize=16)
dev.new(height = 10, width = 14, noRStudioGD = T, dpi=80, pointsize=16)
par(mfrow = c(2,2), mar=c(5,5,1,1), mgp=c(2.5,1,0))

hist(sp_abund$total_abund, breaks=seq(0,500,by=10), main = "", ylab="Number of species", xlab="",las=1, col="grey80")
title(xlab="Abundance", mgp=c(2,1,0))
title(main="(a) All Data", mgp=c(2,1,0), adj = 0, font.main = 1)

hist(sp_abund$total_abund[sp_abund$total_abund<100], breaks=seq(0,100,by=10), main="", ylab="Number of species", xlab="",las=1, col="grey80")
title(xlab="Abundance", mgp=c(2,1,0))
title(main="(b) Less than 100", mgp=c(2,1,0), adj = 0, font.main = 1)

hist(sp_abund$total_abund[sp_abund$total_abund<20], breaks=seq(0,20,by=1), main="", ylab="Number of species", xlab="",las=1, col="grey80", adj = 0)
title(xlab="Abundance", mgp=c(2,1,0))
title(main="(c) Less than 20", mgp=c(2,1,0), adj = 0,font.main = 1)

hist(sp_abund$total_abund[sp_abund$total_abund<10], breaks=seq(0,10,by=1), main="", ylab="Number of species", xlab="",las=1, col="grey80", adj = 0)
title(xlab="Abundance", mgp=c(2,1,0))
title(main="(d) Less than 10", mgp=c(2,1,0), adj = 0, font.main = 1)

less3 <- sp_abund$species[sp_abund$total_abund<3]
less9 <- sp_abund$species[sp_abund$total_abund<9]


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





