

# 02_analyse_data

#use this if starting fresh
#load("04_workspaces/processed_data.RData")
#already analysed data
load("04_workspaces/analysed_data.RData")

library("MASS")
library(AICcmodavg)

# species matrix into data frame for processing
sp_div2 <- as.data.frame(sp_div2)
#save.image("04_workspaces/analysed_data.RData")

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

#save.image("04_workspaces/analysed_data.RData")

#data fully processed and ready to process 9th May 2023 in sum_dat

# species richness as a function of fire category (generalized linear model with negative binomial linear structure). Statistical model = no negative results because count data. 
# Effect of fire can vary depending on location - fire cat*location have no effect on species richness, fire cat AND location have no effect on species richness.
m1_a <- glm.nb(sp_rich~fire_cat*location,data = sum_dat)
m1_b <- glm.nb(sp_rich~fire_cat+location,data = sum_dat)
m1_c <- glm.nb(sp_rich~fire_cat,data = sum_dat)
m1_d <- glm.nb(sp_rich~1,data = sum_dat)

#Likelihood ratio test - compares more complicated model (interaction model) with simpler model 
anova(m1_a, m1_b, test = "F") 
anova(m1_b, m1_c, test = "F") 
anova(m1_c, m1_d, test = "F") 

summary(m1_a);anova(m1_a) # p value for interaction term = 0.89 (fire x location)
summary(m1_b);anova(m1_b) # p value for location term = 0.60
summary(m1_c);anova(m1_c) # p value for fire term = 0.48
AICc(m1_a) #100.3
AICc(m1_b) #85.3
AICc(m1_c) #80.55
AICc(m1_d) #74.7

#m1.set<-list("fire x location"= m1_a, "location"= m1_b, "fire"= m1_c, "null"= m1_d)
#m1.tab<-aictab(cand.set = m1.set, second.ord = T, sort = T)

#write.table(m1.tab,file="m1_tab.txt", quote = F, sep = "\t", row.names=F)

# plot fire+location for abund_25 + shannons

#data frame for fire only predictions

fireonly.pr<-data.frame(fire_cat = factor(levels(sum_dat$fire_cat),levels = levels(sum_dat$fire_cat)))

# fire predictions for models to plot graphs, use m1_c

m1_c.pr<-predict(object = m1_c, newdata = fireonly.pr,type = "response", se.fit = T)
m1_c.pr2<-data.frame(fireonly.pr)
m1_c.pr2$fit<-m1_c.pr$fit
m1_c.pr2$se<-m1_c.pr$se
m1_c.pr2$lci<-m1_c.pr$fit-(m1_c.pr2$se*1.96)
m1_c.pr2$uci<-m1_c.pr$fit+(m1_c.pr2$se*1.96)

#data frame for location only predictions

locationonly.pr<-data.frame(location = factor(levels(sum_dat$location),levels = levels(sum_dat$location)))

# location predictions for models to plot graphs, use m1_b

m1_b.pr<-predict(object = m1_b, newdata = locationonly.pr,type = "response", se.fit = T)
m1_b.pr2<-data.frame(locationonly.pr)
m1_b.pr2$fit<-m1_b.pr$fit
m1_b.pr2$se<-m1_b.pr$se
m1_b.pr2$lci<-m1_b.pr$fit-(m1_b.pr2$se*1.96)
m1_b.pr2$uci<-m1_b.pr$fit+(m1_b.pr2$se*1.96)

#Cum.Wt = cumulative weight as you add each one up. 0.95 = 95% held in _

# all insignificant p values, fire cat and location no effect on sp richness. AIC shows all models with 2 point difference = significant

#save.image("04_workspaces/analysed_data.RData")

head(sum_dat, 6);dim(sum_dat)

summary(sum_dat$simps_ind2) 

# simpson's index as a function of fire category (generalized linear model)
m2_a <- glm(simps_ind2~fire_cat*location,data = sum_dat, family = "Gamma")
m2_b <- glm(simps_ind2~fire_cat+location,data = sum_dat, family = "Gamma")
m2_c <- glm(simps_ind2~fire_cat,data = sum_dat, family = "Gamma")
m2_d <- glm(simps_ind2~1,data = sum_dat, family = "Gamma")

#Likelihood ratio test - compares more complicated model (interaction model) with simpler model 
anova(m2_a, m2_b, test = "F") 
anova(m2_b, m2_c, test = "F") 
anova(m2_c, m2_d, test = "F") 

summary(m2_a);anova(m2_a) # p value for interaction term = 0.55 (fire x location)
summary(m2_b);anova(m2_b) # p value for location term = 0.10
summary(m2_c);anova(m2_c) # p value for fire term = 0.06 # close to significant, plot fire term
AICc(m2_a) #73.9
AICc(m2_b) #60.7
AICc(m2_c) #59.2
AICc(m2_d) #58.5

#m2.set<-list("fire x location"= m2_a, "location"= m2_b, "fire"= m2_c, "null"= m2_d)
#m2.tab<-aictab(cand.set = m2.set, second.ord = T, sort = T)

#write.table(m2.tab,file="m2_tab.txt", quote = F, sep = "\t", row.names=F)

# all insignificant p values but fire term close to significant at p = 0.06 so plot fire term. AIC models no significant diff except between a and b. c is within 2 AICs of null, c holds 35% of weight.
summary(m2_c)

# data frame for fire only predictions
fireonly.pr<-data.frame(fire_cat = factor(levels(sum_dat$fire_cat),levels = levels(sum_dat$fire_cat)))

# fire predictions for models to plot graphs

m2_c.pr<-predict(object = m2_c, newdata = fireonly.pr,type = "response", se.fit = T)
m2_c.pr2<-data.frame(fireonly.pr)
m2_c.pr2$fit<-m2_c.pr$fit
m2_c.pr2$se<-m2_c.pr$se
m2_c.pr2$lci<-m2_c.pr$fit-(m2_c.pr2$se*1.96)
m2_c.pr2$uci<-m2_c.pr$fit+(m2_c.pr2$se*1.96)

#data frame for location only predictions

locationonly.pr<-data.frame(location = factor(levels(sum_dat$location),levels = levels(sum_dat$location)))

# location predictions for models to plot graphs, use m2_b

m2_b.pr<-predict(object = m2_b, newdata = locationonly.pr,type = "response", se.fit = T)
m2_b.pr2<-data.frame(locationonly.pr)
m2_b.pr2$fit<-m2_b.pr$fit
m2_b.pr2$se<-m2_b.pr$se
m2_b.pr2$lci<-m2_b.pr$fit-(m2_b.pr2$se*1.96)
m2_b.pr2$uci<-m2_b.pr$fit+(m2_b.pr2$se*1.96)

# abund_25 as a function of fire category (generalized linear model with negative binomial linear structure)

m3_a <- glm.nb(abund_25~fire_cat*location,data = sum_dat)
m3_b <- glm.nb(abund_25~fire_cat+location,data = sum_dat)
m3_c <- glm.nb(abund_25~fire_cat,data = sum_dat)
m3_d <- glm.nb(abund_25~1,data = sum_dat)

summary(m3_a);anova(m3_a) # p value for interaction term = 0.12
summary(m3_b);anova(m3_b) # p value for location term = 0.0114
summary(m3_c);anova(m3_c) # p value for fire term = 0.075
AICc(m3_a) #62.8
AICc(m3_b) #51.8
AICc(m3_c) #53.1
AICc(m3_d) #50.4

# location signficant and fire term close to significant. Plot fire cat = c

#m3.set<-list("fire x location"= m3_a, "location"= m3_b, "fire"= m3_c, "null"= m3_d)
#m3.tab<-aictab(cand.set = m3.set, second.ord = T, sort = T)

#write.table(m3.tab,file="m3_tab.txt", quote = F, sep = "\t", row.names=F)

#predictions for models to plot graphs
fireonly.pr<-data.frame(fire_cat = factor(levels(sum_dat$fire_cat),levels = levels(sum_dat$fire_cat)))
m3_c.pr<-predict(object = m3_c, newdata = fireonly.pr,type = "response", se.fit = T)
m3_c.pr2<-data.frame(fireonly.pr)
m3_c.pr2$fit<-m3_c.pr$fit
m3_c.pr2$se<-m3_c.pr$se
m3_c.pr2$lci<-m3_c.pr$fit-(m3_c.pr2$se*1.96)
m3_c.pr2$uci<-m3_c.pr$fit+(m3_c.pr2$se*1.96)

#data frame for location only predictions

summary(m3_b)

locationfire.pr<-data.frame(location = rep(factor(levels(sum_dat$location),levels = levels(sum_dat$location)),3),fire_cat = c(rep(levels(sum_dat$fire_cat)[1],2),rep(levels(sum_dat$fire_cat)[2],2),rep(levels(sum_dat$fire_cat)[3],2)))

# location predictions for models to plot graphs, use m3_b

m3_b.pr<-predict(object = m3_b, newdata = locationfire.pr,type = "link", se.fit = T)
m3_b.pr2<-data.frame(locationfire.pr)
m3_b.pr2$fit<-m3_b.pr$fit
m3_b.pr2$se<-m3_b.pr$se
m3_b.pr2$lci<-m3_b.pr$fit-(m3_b.pr2$se*1.96)
m3_b.pr2$uci<-m3_b.pr$fit+(m3_b.pr2$se*1.96)

m3_b.pr2$fit<-exp(m3_b.pr2$fit)
m3_b.pr2$se<-exp(m3_b.pr2$se)
m3_b.pr2$lci<-exp(m3_b.pr2$lci)
m3_b.pr2$uci<-exp(m3_b.pr2$uci)

# abund_5 as a function of fire category (generalized linear model with negative binomial linear structure)
m4_a <- glm.nb(abund_5~fire_cat*location,data = sum_dat)
m4_b <- glm.nb(abund_5~fire_cat+location,data = sum_dat)
m4_c <- glm.nb(abund_5~fire_cat,data = sum_dat)
m4_d <- glm.nb(abund_5~1,data = sum_dat)

summary(m4_a);anova(m4_a) # p value for interaction term = 0.54
summary(m4_b);anova(m4_b) # p value for location term = 0.44
summary(m4_c);anova(m4_c) # p value for fire term = 0.60
AICc(m4_a) #92.4
AICc(m4_b) #78.5
AICc(m4_c) #74.0
AICc(m4_d) #67.7

# all p values insignificant, abund_5 no effect on fire cat. AIC models all significant with each other but none better than null.

#m4.set<-list("fire x location"= m4_a, "location"= m4_b, "fire"= m4_c, "null"= m4_d)
#m4.tab<-aictab(cand.set = m4.set, second.ord = T, sort = T)

#write.table(m4.tab,file="m4_tab.txt", quote = F, sep = "\t", row.names=F)

# abund_5 predictions for models to plot graphs
fireonly.pr<-data.frame(fire_cat = factor(levels(sum_dat$fire_cat),levels = levels(sum_dat$fire_cat)))
m4_c.pr<-predict(object = m4_c, newdata = fireonly.pr,type = "response", se.fit = T)
m4_c.pr2<-data.frame(fireonly.pr)
m4_c.pr2$fit<-m4_c.pr$fit
m4_c.pr2$se<-m4_c.pr$se
m4_c.pr2$lci<-m4_c.pr$fit-(m4_c.pr2$se*1.96)
m4_c.pr2$uci<-m4_c.pr$fit+(m4_c.pr2$se*1.96)

#data frame for location only predictions

locationonly.pr<-data.frame(location = factor(levels(sum_dat$location),levels = levels(sum_dat$location)))

# location predictions for models to plot graphs, use m2_b

m4_b.pr<-predict(object = m4_b, newdata = locationonly.pr,type = "response", se.fit = T)
m4_b.pr2<-data.frame(locationonly.pr)
m4_b.pr2$fit<-m4_b.pr$fit
m4_b.pr2$se<-m4_b.pr$se
m4_b.pr2$lci<-m4_b.pr$fit-(m4_b.pr2$se*1.96)
m4_b.pr2$uci<-m4_b.pr$fit+(m4_b.pr2$se*1.96)

summary(sum_dat$shann_ind)

#species diversity as shannon's index:
# shann_ind as a function of fire category (generalized linear model)
m5_a <- glm(shann_ind~fire_cat*location,data = sum_dat, family = "Gamma")
m5_b <- glm(shann_ind~fire_cat+location,data = sum_dat, family = "Gamma")
m5_c <- glm(shann_ind~fire_cat,data = sum_dat, family = "Gamma")
m5_d <- glm(shann_ind~1,data = sum_dat, family = "Gamma")

#Likelihood ratio test - compares more complicated model (interaction model) with simpler model 
anova(m5_a, m5_b, test = "F") 
anova(m5_b, m5_c, test = "F") 
anova(m5_c, m5_d, test = "F") 

summary(m5_a);anova(m5_a) # p value for interaction term = 0.36 (fire x location)
summary(m5_b);anova(m5_b) # p value for location term = 0.03
summary(m5_c);anova(m5_c) # p value for fire term = 0.02 
AICc(m5_a) #13.5
AICc(m5_b) #1.8
AICc(m5_c) #3.7
AICc(m5_d) #5.3

# location and fire both significant. AIC models all significant with each other, both location and fire better than null. Plot both location and fire?

#m5.set<-list("fire x location"= m5_a, "location"= m5_b, "fire"= m5_c, "null"= m5_d)
#m5.tab<-aictab(cand.set = m5.set, second.ord = T, sort = T)

#write.table(m5.tab,file="m5_tab.txt", quote = F, sep = "\t", row.names=F)

#predictions for models to plot graphs
fireonly.pr<-data.frame(fire_cat = factor(levels(sum_dat$fire_cat),levels = levels(sum_dat$fire_cat)))
m5_c.pr<-predict(object = m5_c, newdata = fireonly.pr,type = "response", se.fit = T)
m5_c.pr2<-data.frame(fireonly.pr)
m5_c.pr2$fit<-m5_c.pr$fit
m5_c.pr2$se<-m5_c.pr$se
m5_c.pr2$lci<-m5_c.pr$fit-(m5_c.pr2$se*1.96)
m5_c.pr2$uci<-m5_c.pr$fit+(m5_c.pr2$se*1.96)

#data frame for location only predictions

summary(m5_b)

locationfire.pr<-data.frame(location = rep(factor(levels(sum_dat$location),levels = levels(sum_dat$location)),3),fire_cat = c(rep(levels(sum_dat$fire_cat)[1],2),rep(levels(sum_dat$fire_cat)[2],2),rep(levels(sum_dat$fire_cat)[3],2)))

# location predictions for models to plot graphs, use m5_b

m5_b.pr<-predict(object = m5_b, newdata = locationfire.pr,type = "link", se.fit = T)
m5_b.pr2<-data.frame(locationfire.pr)
m5_b.pr2$fit<-m5_b.pr$fit
m5_b.pr2$se<-m5_b.pr$se
m5_b.pr2$lci<-m5_b.pr$fit-(m5_b.pr2$se*1.96)
m5_b.pr2$uci<-m5_b.pr$fit+(m5_b.pr2$se*1.96)

#dev.new(width=10, height=5, dpi=80, pointsize=16, noRStudioGD = T)
#par(mfrow=c(1,2), mar=c(4,4,1,1), mgp=c(1.75,0.8,0))

# don't model pres_5 because all 1s, pres_25 not enough 0s

# species diversity as even2:
# even2 as a function of fire category (generalized linear model)
m6_a <- glm(even2~fire_cat*location,data = sum_dat, family = "Gamma")
m6_b <- glm(even2~fire_cat+location,data = sum_dat, family = "Gamma")
m6_c <- glm(even2~fire_cat,data = sum_dat, family = "Gamma")
m6_d <- glm(even2~1,data = sum_dat, family = "Gamma")

#Likelihood ratio test - compares more complicated model (interaction model) with simpler model 
anova(m6_a, m6_b, test = "F") 
anova(m6_b, m6_c, test = "F") 
anova(m6_c, m6_d, test = "F") 

summary(m6_a);anova(m6_a) # p value for interaction term = 0.55 (fire x location)
summary(m6_b);anova(m6_b) # p value for location term = 0.10
summary(m6_c);anova(m6_c) # p value for fire term = 0.06 # close to significant, plot fire term
AICc(m6_a) # -16.7
AICc(m6_b) # -27.4
AICc(m6_c) # -27.7
AICc(m6_d) # -27.4

# all insignificant but fire term close to significant so plot. AIC b, c and d all extremely close, b same as null term so plot?

#predictions for models to plot graphs
fireonly.pr<-data.frame(fire_cat = factor(levels(sum_dat$fire_cat),levels = levels(sum_dat$fire_cat)))
m6_c.pr<-predict(object = m6_c, newdata = fireonly.pr,type = "response", se.fit = T)
m6_c.pr2<-data.frame(fireonly.pr)
m6_c.pr2$fit<-m6_c.pr$fit
m6_c.pr2$se<-m6_c.pr$se
m6_c.pr2$lci<-m6_c.pr$fit-(m6_c.pr2$se*1.96)
m6_c.pr2$uci<-m6_c.pr$fit+(m6_c.pr2$se*1.96)

#m6.set<-list("fire x location"= m6_a, "location"= m6_b, "fire"= m6_c, "null"= m6_d)
#m6.tab<-aictab(cand.set = m6.set, second.ord = T, sort = T)

#m1.tab2<-data.frame(response="Species Richness",m1.tab)
#m2.tab2<-data.frame(response="Simpson's Index",m2.tab)
#m3.tab2<-data.frame(response="Abund_25",m3.tab)
#m4.tab2<-data.frame(response="Abund_5",m4.tab)
#m5.tab2<-data.frame(response="Shannon's Index",m5.tab)
#m6.tab2<-data.frame(response="Evenness",m6.tab)

#combi.tab<-rbind(m1.tab2, m2.tab2, m3.tab2, m4.tab2, m5.tab2, m6.tab2)

#write.table(combi.tab,file="megatable.txt", quote = F, sep = "\t", row.names=F)

# plots for strongest fit

dev.new(width=7, height=9, dpi=80, pointsize=16, noRStudioGD = T)
par(mfrow=c(3,2), mar=c(4.5,4,1,1), mgp=c(2.8,0.8,0), oma=c(0,0,1,6))

# species richness
plot(c(1:3),m1_c.pr2$fit, xlim=c(0.5,3.5), pch=20, xaxt="n",ylim= c((min(m1_c.pr2$lci)),max(m1_c.pr2$uci)),ylab="Species Richness",xlab="", las = 1, cex = 2.5)
arrows(c(1:3),m1_c.pr2$lci,c(1:3),m1_c.pr2$uci,length=0.03,code=3,angle=90)
axis(1,at=c(1:3),labels=F)
axis(1,at=c(0.8,2,3.2),labels=m1_c.pr2$fire_cat,tick=F)
title(mgp=c(2.3,0.8,0),xlab="Fire Category")
m1.tab2
text(2.1,max(m1_c.pr2$uci),as.expression(bquote(Delta~"AICc ="~.(paste(round(m1.tab2$Delta_AICc[m1.tab2$Modnames=="fire"]-m1.tab2$Delta_AICc[m1.tab2$Modnames=="null"],2),sep="")))),adj=0,col="red")
mtext(text="(a)", side = 3, line = 0.5, adj = 0, cex = 1)

 # simps diversity index, no changes required
plot(c(1:3),m2_c.pr2$fit, xlim=c(0.5,3.5), pch=20, xaxt="n",ylim= c((min(m2_c.pr2$lci)),max(m2_c.pr2$uci)),ylab="Simpson's Diversity Index",xlab="", las = 1, cex = 2.5)
arrows(c(1:3),m2_c.pr2$lci,c(1:3),m2_c.pr2$uci,length=0.03,code=3,angle=90)
axis(1,at=c(1:3),labels=F)
axis(1,at=c(0.8,2,3.2),labels=m2_c.pr2$fire_cat,tick=F)
title(mgp=c(2.3,0.8,0),xlab="Fire Category")
text(2.1,max(m2_c.pr2$uci),as.expression(bquote(Delta~"AICc ="~.(paste(round(m2.tab2$Delta_AICc[m2.tab2$Modnames=="fire"]-m2.tab2$Delta_AICc[m2.tab2$Modnames=="null"],2),sep="")))),adj=0,col="red")
mtext(text="(b)", side = 3, line = 0.5, adj = 0, cex = 1)

par(xpd=NA)
legend(x=4,y=max(m2_c.pr2$uci)+0.1, title = "Sites", legend = c("Fire only", "Hincks","Pinks"), pt.cex = 1.5, pch = c(16, 15, 17), bty = "n", title.adj=0)
par(xpd=F)

# shann_ind plots for fire+location

plot(c(1:3)-0.35,m5_b.pr2$fit[m5_b.pr2$location=="Hincks"], xlim=c(0.5,3.5), pch=15, xaxt="n",ylim= c((min(m5_b.pr2$lci)),max(m5_b.pr2$uci)+0.03),ylab="Shannon's Index",xlab="", las = 1, cex = 1.5)
points(c(1:3),m5_b.pr2$fit[m5_b.pr2$location=="Pinks"], xlim=c(0.5,3.5), pch=17, cex = 1.5)
arrows(c(1:3)-0.35,m5_b.pr2$lci[m5_b.pr2$location=="Hincks"],c(1:3)-0.35,m5_b.pr2$uci[m5_b.pr2$location=="Hincks"],length=0.03,code=3,angle=90)
arrows(c(1:3),m5_b.pr2$lci[m5_b.pr2$location=="Pinks"],c(1:3),m5_b.pr2$uci[m5_b.pr2$location=="Pinks"],length=0.03,code=3,angle=90)
axis(1,at=c(1:3)-0.25,labels=F)
axis(1,at=c(0.8,2,3.2),labels=m5_b.pr2$fire_cat[m5_b.pr2$location=="Hincks"],tick=F)
title(mgp=c(2.3,0.8,0),xlab="Fire Category")
text(2,max(m5_b.pr2$uci)+0.03,as.expression(bquote(Delta~"AICc ="~.(paste(round(m5.tab2$Delta_AICc[m5.tab2$Modnames=="fire"]-m5.tab2$Delta_AICc[m5.tab2$Modnames=="null"],2),sep="")))),adj=0,col="dark green")
mtext(text="(c)", side = 3, line = 0.5, adj = 0, cex = 1)

# evenness plots for fire
plot(c(1:3),m6_c.pr2$fit, xlim=c(0.5,3.5), pch=20, xaxt="n",ylim= c((min(m6_c.pr2$lci)),max(m6_c.pr2$uci)),ylab="Evenness",xlab="", las = 1, cex = 2.5)
arrows(c(1:3),m6_c.pr2$lci,c(1:3),m6_c.pr2$uci,length=0.03,code=3,angle=90)
axis(1,at=c(1:3),labels=F)
axis(1,at=c(0.8,2,3.2),labels=m6_c.pr2$fire_cat,tick=F)
title(mgp=c(2.3,0.8,0),xlab="Fire category")
text(2,max(m6_c.pr2$uci),as.expression(bquote(Delta~"AICc ="~.(paste(round(m6.tab2$Delta_AICc[m6.tab2$Modnames=="fire"]-m6.tab2$Delta_AICc[m6.tab2$Modnames=="null"],2),sep="")))),adj=0,col="dark green")
mtext(text="(d)", side = 3, line = 0.5, adj = 0, cex = 1)

# abund_5 plot
plot(c(1:3),m4_c.pr2$fit, xlim=c(0.5,3.5), pch=20, xaxt="n",ylim= c((min(m4_c.pr2$lci)),max(m4_c.pr2$uci)),ylab="5% of max abundance",xlab="", las = 1, cex = 2.5)
arrows(c(1:3),m4_c.pr2$lci,c(1:3),m4_c.pr2$uci,length=0.03,code=3,angle=90)
axis(1,at=c(1:3),labels=F)
axis(1,at=c(0.8,2,3.2),labels=m4_c.pr2$fire_cat,tick=F)
title(mgp=c(2.3,0.8,0),xlab="Fire Category")
text(2.1,max(m4_c.pr2$uci),as.expression(bquote(Delta~"AICc ="~.(paste(round(m4.tab2$Delta_AICc[m4.tab2$Modnames=="fire"]-m4.tab2$Delta_AICc[m4.tab2$Modnames=="null"],2),sep="")))),adj=0,col="red")
mtext(text="(e)", side = 3, line = 0.5, adj = 0, cex = 1)

# abund_25 plot fire*location

plot(c(1:3)-0.35,m3_b.pr2$fit[m3_b.pr2$location=="Hincks"], xlim=c(0.5,3.5), pch=15, xaxt="n",ylim= c((min(m3_b.pr2$lci)),max(m3_b.pr2$uci)),ylab="Lowest 25%",xlab="", las = 1, cex = 1.5)
points(c(1:3),m3_b.pr2$fit[m3_b.pr2$location=="Pinks"], xlim=c(0.5,3.5), pch=17, cex = 1.5)
arrows(c(1:3)-0.35,m3_b.pr2$lci[m3_b.pr2$location=="Hincks"],c(1:3)-0.35,m3_b.pr2$uci[m3_b.pr2$location=="Hincks"],length=0.03,code=3,angle=90)
arrows(c(1:3),m3_b.pr2$lci[m3_b.pr2$location=="Pinks"],c(1:3),m3_b.pr2$uci[m3_b.pr2$location=="Pinks"],length=0.03,code=3,angle=90)
axis(1,at=c(1:3)-0.25,labels=F)
axis(1,at=c(0.8,2,3.2),labels=m3_b.pr2$fire_cat[m3_b.pr2$location=="Hincks"],tick=F)
title(mgp=c(2.3,0.8,0),xlab="Fire Category")
text(2.1,max(m3_b.pr2$uci),as.expression(bquote(Delta~"AICc ="~.(paste(round(m3.tab2$Delta_AICc[m3.tab2$Modnames=="fire"]-m3.tab2$Delta_AICc[m3.tab2$Modnames=="null"],2),sep="")))),adj=0,col="red")
mtext(text="(f)", side = 3, line = 0.5, adj = 0, cex = 1)


### start of casual graph plots


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




