# ------------------------------------ #
# ---------- RARE REPTILES  ---------- #
# ------------------------------------ #

### Analysis fire effects on rarity and dominance in a reptile community 

# ---------- 02_analyse_data  ---------- #

### Script authors: Amber Lim & Annabel Smith 

# use this if starting fresh
# load("04_workspaces/processed_data.RData")

# already analysed data
load("04_workspaces/analysed_data.RData")
# save.image("04_workspaces/analysed_data.RData")

library("MASS")
library(AICcmodavg)
library(mgcv)

# Load functions:
invisible(lapply(paste("02_functions/",dir("02_functions"),sep=""), function(x) source(x)))

# data fully processed and ready to analyse 9th May 2023 in sum_dat

head(sum_dat, 3);dim(sum_dat)

# Data distributions and model types: ----
# whole community:
summary(sum_dat$sp_rich) # neg bin
summary(sum_dat$simps_ind2) # gamma
summary(sum_dat$shann_ind) # gamma 
summary(sum_dat$even2) # beta

summary(sum_dat$bp_ind) # beta
summary(sum_dat$fa.all) # gamma

# rare species:
summary(sum_dat$abund_25) # neg bin
summary(sum_dat$abund_5) # neg bin
summary(sum_dat$sr_25) # neg bin
summary(sum_dat$sr_5) # neg bin

# Data distributions:
dev.new(width=6,height=12, dpi=70, pointsize=16,noRStudioGD = T)
par(mfrow=c(5,2),mar=c(4,4,1,1), mgp=c(2.5,1,0))

# whole community:
hist(sum_dat$sp_rich, xlab="",main="species richness", font.main=1, las=1) 
hist(sum_dat$simps_ind2, xlab="",main="simpson's index", font.main=1, las=1) 
hist(sum_dat$shann_ind, xlab="",main="shannon's index", font.main=1, las=1)
hist(sum_dat$even2, xlim=c(0,1), xlab="",main="evenness", font.main=1, las=1) 

hist(sum_dat$bp_ind, xlim=c(0,1), xlab="",main="berger-parker", font.main=1, las=1) 
hist(sum_dat$fa.all, xlab="",main="fisher's alpha", font.main=1, las=1)

# rare species:
hist(sum_dat$abund_25, xlab="",main="abundance 25", font.main=1, las=1) 
hist(sum_dat$abund_5, xlab="",main="abundance 5", font.main=1, las=1) 
hist(sum_dat$sr_25, xlab="",main="richness 25", font.main=1, las=1) 
hist(sum_dat$sr_5, xlab="",main="richness 5", font.main=1, las=1) 

# ----

#### MODELS

# ** species richness (negative binomial GLM) ----

m1_a <- glm.nb(sp_rich~fire_cat*location,data = sum_dat)
m1_b <- glm.nb(sp_rich~fire_cat+location,data = sum_dat)
m1_c <- glm.nb(sp_rich~fire_cat,data = sum_dat)
m1_d <- glm.nb(sp_rich~1,data = sum_dat)

#Likelihood ratio test
anova(m1_a, m1_b, test = "F") 
anova(m1_b, m1_c, test = "F") 
anova(m1_c, m1_d, test = "F") 

summary(m1_a);anova(m1_a) 
summary(m1_b);anova(m1_b) 
summary(m1_c);anova(m1_c) 

AICc(m1_a) #100.3
AICc(m1_b) #85.3
AICc(m1_c) #80.55
AICc(m1_d) #74.7

m1.set<-list("fire x location"= m1_a, "location"= m1_b, "fire"= m1_c, "null"= m1_d)
m1.tab<-aictab(cand.set = m1.set, second.ord = T, sort = T)

# fire only predictions (m1_c):

# For neg bin models, predict on the link scale and backtransform after, to ensure CIs are bounded by zero:

fireonly.pr<-data.frame(fire_cat = factor(levels(sum_dat$fire_cat),levels = levels(sum_dat$fire_cat)))

m1_c.pr<-predict(object = m1_c, newdata = fireonly.pr,type = "link", se.fit = T)
m1_c.pr2<-data.frame(fireonly.pr)
m1_c.pr2$fit<-m1_c.pr$fit
m1_c.pr2$se<-m1_c.pr$se
m1_c.pr2$lci<-m1_c.pr$fit-(m1_c.pr2$se*1.96)
m1_c.pr2$uci<-m1_c.pr$fit+(m1_c.pr2$se*1.96)

m1_c.pr2$lci85<-m1_c.pr$fit-(m1_c.pr2$se*1.44)
m1_c.pr2$uci85<-m1_c.pr$fit+(m1_c.pr2$se*1.44)

m1_c.pr2$fit<-exp(m1_c.pr2$fit)
m1_c.pr2$se<-exp(m1_c.pr2$se)
m1_c.pr2$lci<-exp(m1_c.pr2$lci)
m1_c.pr2$uci<-exp(m1_c.pr2$uci)

m1_c.pr2$lci85<-exp(m1_c.pr2$lci85)
m1_c.pr2$uci85<-exp(m1_c.pr2$uci85)

# save.image("04_workspaces/analysed_data.RData")

head(sum_dat, 6);dim(sum_dat)

# ---- 

# ** simpson's index (Gamma GLM) ---- 

summary(sum_dat$simps_ind2) 

m2_a <- glm(simps_ind2~fire_cat*location,data = sum_dat, family = "Gamma")
m2_b <- glm(simps_ind2~fire_cat+location,data = sum_dat, family = "Gamma")
m2_c <- glm(simps_ind2~fire_cat,data = sum_dat, family = "Gamma")
m2_d <- glm(simps_ind2~1,data = sum_dat, family = "Gamma")

# Likelihood ratio test 
anova(m2_a, m2_b, test = "F") 
anova(m2_b, m2_c, test = "F") 
anova(m2_c, m2_d, test = "F") 

summary(m2_a);anova(m2_a) 
summary(m2_b);anova(m2_b) 
summary(m2_c);anova(m2_c) 

AICc(m2_a) #73.9
AICc(m2_b) #60.7
AICc(m2_c) #59.2
AICc(m2_d) #58.5

# simpson's index model set:
m2.set<-list("fire x location"= m2_a, "location"= m2_b, "fire"= m2_c, "null"= m2_d)
m2.tab<-aictab(cand.set = m2.set, second.ord = T, sort = T)

# fire only predictions (m2_c)

summary(m2_c)
fireonly.pr

m2_c.pr<-predict(object = m2_c, newdata = fireonly.pr,type = "response", se.fit = T)
m2_c.pr2<-data.frame(fireonly.pr)
m2_c.pr2$fit<-m2_c.pr$fit
m2_c.pr2$se<-m2_c.pr$se
m2_c.pr2$lci<-m2_c.pr$fit-(m2_c.pr2$se*1.96)
m2_c.pr2$uci<-m2_c.pr$fit+(m2_c.pr2$se*1.96)

m2_c.pr2$lci85<-m2_c.pr2$fit-(m2_c.pr2$se*1.44)
m2_c.pr2$uci85<-m2_c.pr2$fit+(m2_c.pr2$se*1.44)

# ---- 

# ** shannon's index (Gamma GLM) ---- 

m5_a <- glm(shann_ind~fire_cat*location,data = sum_dat, family = "Gamma")
m5_b <- glm(shann_ind~fire_cat+location,data = sum_dat, family = "Gamma")
m5_c <- glm(shann_ind~fire_cat,data = sum_dat, family = "Gamma")
m5_d <- glm(shann_ind~1,data = sum_dat, family = "Gamma")

# Likelihood ratio test
anova(m5_a, m5_b, test = "F") 
anova(m5_b, m5_c, test = "F") 
anova(m5_c, m5_d, test = "F") 

summary(m5_a);anova(m5_a) 
summary(m5_b);anova(m5_b) 
summary(m5_c);anova(m5_c) 

AICc(m5_a) #13.5
AICc(m5_b) #1.8
AICc(m5_c) #3.7
AICc(m5_d) #5.3

# shannon's model set:
m5.set<-list("fire x location"= m5_a, "location"= m5_b, "fire"= m5_c, "null"= m5_d)
m5.tab<-aictab(cand.set = m5.set, second.ord = T, sort = T)

# fire + location predictions, m5_b

locationfire.pr<-data.frame(location = rep(factor(levels(sum_dat$location),levels = levels(sum_dat$location)),3),fire_cat = c(rep(levels(sum_dat$fire_cat)[1],2),rep(levels(sum_dat$fire_cat)[2],2),rep(levels(sum_dat$fire_cat)[3],2)))

summary(m5_b)

m5_b.pr<-predict(object = m5_b, newdata = locationfire.pr,type = "response", se.fit = T)
m5_b.pr2<-data.frame(locationfire.pr)
m5_b.pr2$fit<-m5_b.pr$fit
m5_b.pr2$se<-m5_b.pr$se
m5_b.pr2$lci<-m5_b.pr$fit-(m5_b.pr2$se*1.96)
m5_b.pr2$uci<-m5_b.pr$fit+(m5_b.pr2$se*1.96)

m5_b.pr2$lci85<-m5_b.pr2$fit-(m5_b.pr2$se*1.44)
m5_b.pr2$uci85<-m5_b.pr2$fit+(m5_b.pr2$se*1.44)

# ---- 

# ** species evenness (even2) (beta regression): ---- 

# update 14 Jan 2025:

# Should evenness and Berger Parker be modelled with a beta gam, rather than Gamma glm, because the values range between zero and 1? Probably yes, so re-doing this model set. 

# From Ella's thesis: "...GAMM with a beta distribution using the ‘gamm’ function in mgcv (Wood 2017)"

# Must specify ML since we're estimating AICc (the default is fitting by REML)

# even2 (beta regression)
m6_a_beta <- gam(even2~fire_cat*location,data = sum_dat, family = "betar",method="ML")
m6_b_beta <- gam(even2~fire_cat+location,data = sum_dat, family = "betar",method="ML")
m6_c_beta <- gam(even2~fire_cat,data = sum_dat, family = "betar",method="ML")
m6_d_beta <- gam(even2~1,data = sum_dat, family = "betar",method="ML")

summary(m6_a_beta);anova(m6_a_beta) 
summary(m6_b_beta);anova(m6_b_beta) 
summary(m6_c_beta);anova(m6_c_beta) 

# Likelihood ratio test 
anova(m6_a_beta, m6_b_beta, test = "F") 
anova(m6_b_beta, m6_c_beta, test = "F") 
anova(m6_c_beta, m6_d_beta, test = "F") 

AICc(m6_a_beta) # 
AICc(m6_b_beta) # 
AICc(m6_c_beta) # 
AICc(m6_d_beta) # 

# evenness model set
m6.set<-list("fire x location"= m6_a_beta, "location"= m6_b_beta, "fire"= m6_c_beta, "null"= m6_d_beta)

# cannot use aictab for beta models, so need to construct our own:
m6.tab<-data.frame(Modnames=c("fire x location", "location", "fire", "null"))

# K:
m6.tab$K<-c(attributes(logLik.gam(m6_a_beta))$df,attributes(logLik.gam(m6_b_beta))$df,attributes(logLik.gam(m6_c_beta))$df,attributes(logLik.gam(m6_d_beta))$df)

# AICc:
m6.tab$AICc<-c(AICc(m6_a_beta), AICc(m6_b_beta), AICc(m6_c_beta), AICc(m6_d_beta))

# Log-likelihood
m6.tab$LL<-round(c(logLik.gam(m6_a_beta),
                   logLik.gam(m6_b_beta),
                   logLik.gam(m6_c_beta),
                   logLik.gam(m6_d_beta)),2)

# Order:
m6.tab<-m6.tab[order(m6.tab$AICc),]
rownames(m6.tab)<-1:nrow(m6.tab)

# Delta AICc

m6.tab$Delta_AICc<-round(c(0,m6.tab$AICc[-1] - m6.tab$AICc[1]),2)

# Model likelihood: the relative likelihood of the model given the data (exp(-0.5*delta[i]))

m6.tab$ModelLik<-round(exp(-0.5*m6.tab$Delta_AICc),2)

# AICc weight

# https://brianomeara.info/aic.html: Akaike weights are the relative likelihood divided by the sum of these values across all models.

m6.tab$AICcWt<-round(m6.tab$ModelLik/sum(m6.tab$ModelLik),2)

# Cumulative weight

m6.tab$Cum.Wt<-round(cumsum(m6.tab$AICcWt),2)

# Reorder so they align with the others (note LL needs to be added before the change in order):
m6.tab<-m6.tab[,colnames(data.frame(m1.tab))]

# save.image("04_workspaces/analysed_data.RData")

# fire + location predictions, m6_b_beta

summary(m6_b_beta)

m6_b.pr<-predict(object = m6_b_beta, newdata = locationfire.pr,type = "response", se.fit = T)
m6_b.pr2<-data.frame(locationfire.pr)
m6_b.pr2$fit<-m6_b.pr$fit
m6_b.pr2$se<-m6_b.pr$se
m6_b.pr2$lci<-m6_b.pr$fit-(m6_b.pr2$se*1.96)
m6_b.pr2$uci<-m6_b.pr$fit+(m6_b.pr2$se*1.96)

m6_b.pr2$lci85<-m6_b.pr2$fit-(m6_b.pr2$se*1.44)
m6_b.pr2$uci85<-m6_b.pr2$fit+(m6_b.pr2$se*1.44)

# ----

# abundance of RARE SPECIES (negative binomial GLM)

# ** abund_25 (negative binomial GLM) ---- 

head(sum_dat, 2);dim(sum_dat)

m3_a <- glm.nb(abund_25~fire_cat*location,data = sum_dat)
m3_b <- glm.nb(abund_25~fire_cat+location,data = sum_dat)
m3_c <- glm.nb(abund_25~fire_cat,data = sum_dat)
m3_d <- glm.nb(abund_25~1,data = sum_dat)

summary(m3_a);anova(m3_a) 
summary(m3_b);anova(m3_b) 
summary(m3_c);anova(m3_c) 
AICc(m3_a) #62.8
AICc(m3_b) #51.8
AICc(m3_c) #53.1
AICc(m3_d) #50.4

# Abund 25 model set:
m3.set<-list("fire x location"= m3_a, "location"= m3_b, "fire"= m3_c, "null"= m3_d)
m3.tab<-aictab(cand.set = m3.set, second.ord = T, sort = T)

# Abund 25 fire + location predictions

summary(m3_b)

# fire + location predictions, m3_b
# For neg bin models, predict on the link scale and backtransform after, to ensure CIs are bounded by zero:

m3_b.pr<-predict(object = m3_b, newdata = locationfire.pr,type = "link", se.fit = T)
m3_b.pr2<-data.frame(locationfire.pr)
m3_b.pr2$fit<-m3_b.pr$fit
m3_b.pr2$se<-m3_b.pr$se
m3_b.pr2$lci<-m3_b.pr$fit-(m3_b.pr2$se*1.96)
m3_b.pr2$uci<-m3_b.pr$fit+(m3_b.pr2$se*1.96)

m3_b.pr2$lci85<-m3_b.pr2$fit-(m3_b.pr2$se*1.44)
m3_b.pr2$uci85<-m3_b.pr2$fit+(m3_b.pr2$se*1.44)

m3_b.pr2$fit<-exp(m3_b.pr2$fit)
m3_b.pr2$se<-exp(m3_b.pr2$se)
m3_b.pr2$lci<-exp(m3_b.pr2$lci)
m3_b.pr2$uci<-exp(m3_b.pr2$uci)

m3_b.pr2$lci85<-exp(m3_b.pr2$lci85)
m3_b.pr2$uci85<-exp(m3_b.pr2$uci85)

# ---- 

# ** abund_5 (negative binomial GLM) ---- 

m4_a <- glm.nb(abund_5~fire_cat*location,data = sum_dat)
m4_b <- glm.nb(abund_5~fire_cat+location,data = sum_dat)
m4_c <- glm.nb(abund_5~fire_cat,data = sum_dat)
m4_d <- glm.nb(abund_5~1,data = sum_dat)

summary(m4_a);anova(m4_a) 
summary(m4_b);anova(m4_b) 
summary(m4_c);anova(m4_c) 

AICc(m4_a) #92.4
AICc(m4_b) #78.5
AICc(m4_c) #74.0
AICc(m4_d) #67.7

# abund 5 model set:
m4.set<-list("fire x location"= m4_a, "location"= m4_b, "fire"= m4_c, "null"= m4_d)
m4.tab<-aictab(cand.set = m4.set, second.ord = T, sort = T)

# abund_5 predictions 
# For neg bin models, predict on the link scale and backtransform after, to ensure CIs are bounded by zero:
m4_c.pr<-predict(object = m4_c, newdata = fireonly.pr,type = "link", se.fit = T)
m4_c.pr2<-data.frame(fireonly.pr)
m4_c.pr2$fit<-m4_c.pr$fit
m4_c.pr2$se<-m4_c.pr$se
m4_c.pr2$lci<-m4_c.pr$fit-(m4_c.pr2$se*1.96)
m4_c.pr2$uci<-m4_c.pr$fit+(m4_c.pr2$se*1.96)

m4_c.pr2$lci85<-m4_c.pr2$fit-(m4_c.pr2$se*1.44)
m4_c.pr2$uci85<-m4_c.pr2$fit+(m4_c.pr2$se*1.44)

m4_c.pr2$fit<-exp(m4_c.pr2$fit)
m4_c.pr2$se<-exp(m4_c.pr2$se)
m4_c.pr2$lci<-exp(m4_c.pr2$lci)
m4_c.pr2$uci<-exp(m4_c.pr2$uci)

m4_c.pr2$lci85<-exp(m4_c.pr2$lci85)
m4_c.pr2$uci85<-exp(m4_c.pr2$uci85)

# ---- 

# species richness of RARE SPECIES (negative binomial GLM)

head(sum_dat,3); dim(sum_dat)

# Rare species richness, lowest 25 % ----

m7_a <- glm.nb(sr_25~fire_cat*location,data = sum_dat)
m7_b <- glm.nb(sr_25~fire_cat+location,data = sum_dat)
m7_c <- glm.nb(sr_25~fire_cat,data = sum_dat)
m7_d <- glm.nb(sr_25~1,data = sum_dat)

#Likelihood ratio test - compare more complicated model (interaction model) with simpler model 
anova(m7_a, m7_b, test = "F") 
anova(m7_b, m7_c, test = "F") 
anova(m7_c, m7_d, test = "F") 

summary(m7_a);anova(m7_a) # p value for interaction term 
summary(m7_b);anova(m7_b) # p value for location term 
summary(m7_c);anova(m7_c) # p value for fire term 

AICc(m7_a) 
AICc(m7_b)
AICc(m7_c)
AICc(m7_d)

m7.set<-list("fire x location"= m7_a, "location"= m7_b, "fire"= m7_c, "null"= m7_d)
m7.tab<-aictab(cand.set = m7.set, second.ord = T, sort = T)

#data frame for fire only predictions
fireonly.pr

# For neg bin models, predict on the link scale and backtransform after, to ensure CIs are bounded by zero:
m7_c.pr<-predict(object = m7_c, newdata = fireonly.pr,type = "link", se.fit = T)
m7_c.pr2<-data.frame(fireonly.pr)
m7_c.pr2$fit<-m7_c.pr$fit
m7_c.pr2$se<-m7_c.pr$se
m7_c.pr2$lci<-m7_c.pr$fit-(m7_c.pr2$se*1.96)
m7_c.pr2$uci<-m7_c.pr$fit+(m7_c.pr2$se*1.96)

m7_c.pr2$lci85<-m7_c.pr2$fit-(m7_c.pr2$se*1.44)
m7_c.pr2$uci85<-m7_c.pr2$fit+(m7_c.pr2$se*1.44)

m7_c.pr2$fit<-exp(m7_c.pr2$fit)
m7_c.pr2$se<-exp(m7_c.pr2$se)
m7_c.pr2$lci<-exp(m7_c.pr2$lci)
m7_c.pr2$uci<-exp(m7_c.pr2$uci)

m7_c.pr2$lci85<-exp(m7_c.pr2$lci85)
m7_c.pr2$uci85<-exp(m7_c.pr2$uci85)

# ----

# Rare species richness, max 5 % ----

m8_a <- glm.nb(sr_5~fire_cat*location,data = sum_dat)
m8_b <- glm.nb(sr_5~fire_cat+location,data = sum_dat)
m8_c <- glm.nb(sr_5~fire_cat,data = sum_dat)
m8_d <- glm.nb(sr_5~1,data = sum_dat)

#Likelihood ratio test
anova(m8_a, m8_b, test = "F") 
anova(m8_b, m8_c, test = "F") 
anova(m8_c, m8_d, test = "F") 

summary(m8_a);anova(m8_a) # p value for interaction term 
summary(m8_b);anova(m8_b) # p value for location term 
summary(m8_c);anova(m8_c) # p value for fire term 

AICc(m8_a) 
AICc(m8_b)
AICc(m8_c)
AICc(m8_d)

m8.set<-list("fire x location"= m8_a, "location"= m8_b, "fire"= m8_c, "null"= m8_d)
m8.tab<-aictab(cand.set = m8.set, second.ord = T, sort = T)

#data frame for fire only predictions
# For neg bin models, predict on the link scale and backtransform after, to ensure CIs are bounded by zero:

fireonly.pr

m8_c.pr<-predict(object = m8_c, newdata = fireonly.pr,type = "link", se.fit = T)
m8_c.pr2<-data.frame(fireonly.pr)
m8_c.pr2$fit<-m8_c.pr$fit
m8_c.pr2$se<-m8_c.pr$se
m8_c.pr2$lci<-m8_c.pr$fit-(m8_c.pr2$se*1.96)
m8_c.pr2$uci<-m8_c.pr$fit+(m8_c.pr2$se*1.96)

m8_c.pr2$lci85<-m8_c.pr2$fit-(m8_c.pr2$se*1.44)
m8_c.pr2$uci85<-m8_c.pr2$fit+(m8_c.pr2$se*1.44)

m8_c.pr2$fit<-exp(m8_c.pr2$fit)
m8_c.pr2$se<-exp(m8_c.pr2$se)
m8_c.pr2$lci<-exp(m8_c.pr2$lci)
m8_c.pr2$uci<-exp(m8_c.pr2$uci)

m8_c.pr2$lci85<-exp(m8_c.pr2$lci85)
m8_c.pr2$uci85<-exp(m8_c.pr2$uci85)

# ----

# save.image("04_workspaces/analysed_data.RData")

# ---*** MODEL BERGER-PARKER & FISHER'S ALPHA ***--- #

# Whole community only

# Berger-Parker (beta regression) ----
# Berger-Parker ranges from 0 to 1, with 0 indicating complete evenness or equal abundance among all species in the community; model as beta regression

head(sum_dat,3); dim(sum_dat)

# bp (beta regression)
m9_a <- gam(bp_ind~fire_cat*location,data = sum_dat, family = "betar",method="ML")
m9_b <- gam(bp_ind~fire_cat+location,data = sum_dat, family = "betar",method="ML")
m9_c <- gam(bp_ind~fire_cat,data = sum_dat, family = "betar",method="ML")
m9_d <- gam(bp_ind~1,data = sum_dat, family = "betar",method="ML")

summary(m9_a);anova(m9_a) 
summary(m9_b);anova(m9_b) 
summary(m9_c);anova(m9_c) 

# Likelihood ratio test 
anova(m9_a, m9_b, test = "F") 
anova(m9_b, m9_c, test = "F") 
anova(m9_c, m9_d, test = "F") 

AICc(m9_a) # 
AICc(m9_b) # 
AICc(m9_c) # 
AICc(m9_d) # 

# bp model set

# cannot use aictab for beta models, so need to construct our own:
m9.tab<-data.frame(Modnames=c("fire x location", "location", "fire", "null"))

# K:
m9.tab$K<-c(attributes(logLik.gam(m9_a))$df,attributes(logLik.gam(m9_b))$df,attributes(logLik.gam(m9_c))$df,attributes(logLik.gam(m9_d))$df)

# AICc:
m9.tab$AICc<-c(AICc(m9_a), AICc(m9_b), AICc(m9_c), AICc(m9_d))

# Log-likelihood
m9.tab$LL<-round(c(logLik.gam(m9_a),
                   logLik.gam(m9_b),
                   logLik.gam(m9_c),
                   logLik.gam(m9_d)),2)

# Order:
m9.tab<-m9.tab[order(m9.tab$AICc),]
rownames(m9.tab)<-1:nrow(m9.tab)

# Delta AICc

m9.tab$Delta_AICc<-round(c(0,m9.tab$AICc[-1] - m9.tab$AICc[1]),2)

# Model likelihood: the relative likelihood of the model given the data (exp(-0.5*delta[i]))

m9.tab$ModelLik<-round(exp(-0.5*m9.tab$Delta_AICc),2)

# AICc weight

# https://brianomeara.info/aic.html: Akaike weights are the relative likelihood divided by the sum of these values across all models.

m9.tab$AICcWt<-round(m9.tab$ModelLik/sum(m9.tab$ModelLik),2)

# Cumulative weight

m9.tab$Cum.Wt<-round(cumsum(m9.tab$AICcWt),2)

# Reorder so they align with the others (note LL needs to be added before the change in order):
m9.tab<-m9.tab[,colnames(data.frame(m1.tab))]

# save.image("04_workspaces/analysed_data.RData")

# fire-only predictions 
m9_c.pr<-predict(object = m9_c, newdata = fireonly.pr,type = "response", se.fit = T)
m9_c.pr2<-data.frame(fireonly.pr)
m9_c.pr2$fit<-m9_c.pr$fit
m9_c.pr2$se<-m9_c.pr$se
m9_c.pr2$lci<-m9_c.pr2$fit-(m9_c.pr2$se*1.96)
m9_c.pr2$uci<-m9_c.pr2$fit+(m9_c.pr2$se*1.96)

m9_c.pr2$lci85<-m9_c.pr2$fit-(m9_c.pr2$se*1.44)
m9_c.pr2$uci85<-m9_c.pr2$fit+(m9_c.pr2$se*1.44)

# ----

# Fisher's alpha (gamma) ----

# Fisher's alpha is positive continuous; model as gamma

head(sum_dat, 3); dim(sum_dat)

m10_a <- glm(fa.all~fire_cat*location,data = sum_dat, family = "Gamma")
m10_b <- glm(fa.all~fire_cat+location,data = sum_dat, family = "Gamma")
m10_c <- glm(fa.all~fire_cat,data = sum_dat, family = "Gamma")
m10_d <- glm(fa.all~1,data = sum_dat, family = "Gamma")

# Likelihood ratio test 
anova(m10_a, m10_b, test = "F") 
anova(m10_b, m10_c, test = "F") 
anova(m10_c, m10_d, test = "F") 

summary(m10_a);anova(m10_a) 
summary(m10_b);anova(m10_b) 
summary(m10_c);anova(m10_c) 

AICc(m10_a) 
AICc(m10_b) 
AICc(m10_c) 
AICc(m10_d)

# fa model set:
m10.set<-list("fire x location"= m10_a, "location"= m10_b, "fire"= m10_c, "null"= m10_d)
m10.tab<-aictab(cand.set = m10.set, second.ord = T, sort = T)

# fire + location predictions, m10_b

summary(m10_b)

m10_b.pr<-predict(object = m10_b, newdata = locationfire.pr,type = "response", se.fit = T)
m10_b.pr2<-data.frame(locationfire.pr)
m10_b.pr2$fit<-m10_b.pr$fit
m10_b.pr2$se<-m10_b.pr$se
m10_b.pr2$lci<-m10_b.pr2$fit-(m10_b.pr2$se*1.96)
m10_b.pr2$uci<-m10_b.pr2$fit+(m10_b.pr2$se*1.96)

m10_b.pr2$lci85<-m10_b.pr2$fit-(m10_b.pr2$se*1.44)
m10_b.pr2$uci85<-m10_b.pr2$fit+(m10_b.pr2$se*1.44)

# ----

# save.image("04_workspaces/analysed_data.RData")

# AICc table: all response variables ----

m1.tab2<-data.frame(response="Species Richness",m1.tab)
m2.tab2<-data.frame(response="Simpson's Index",m2.tab)
m3.tab2<-data.frame(response="Abund_25",m3.tab)
m4.tab2<-data.frame(response="Abund_5",m4.tab)
m5.tab2<-data.frame(response="Shannon's Index",m5.tab)
m6.tab2<-data.frame(response="Evenness",m6.tab)
m7.tab2<-data.frame(response="Richness Lowest 25%",m7.tab)
m8.tab2<-data.frame(response="Richness Max. 5%",m8.tab)

m9.tab2<-data.frame(response="Berger-Parker",m9.tab)
m10.tab2<-data.frame(response="Fisher's Alpha",m10.tab)

# Added m7 and m8 to combi.tab (richness of rare species):
combi.tab<-rbind(m1.tab2, m2.tab2, m3.tab2, m4.tab2, m5.tab2, m6.tab2, m7.tab2, m8.tab2, m9.tab2, m10.tab2)

# write.table(combi.tab,file="AIC_all.txt", quote = F, sep = "\t", row.names=F)

# save.image("04_workspaces/analysed_data.RData")

# ----

#### CONTRASTS (not used in new MS) 

# Set up contrast matrices ----

# fire only contrasts

# There are 3 contrasts for three fire categories (u:m, u:b, m:b)

summary(m1_c) # glm.nb model

fireonly_contrast<-data.frame(contrast=paste(combn(fireonly.pr$fire_cat,2)[1,],combn(fireonly.pr$fire_cat,2)[2,],sep=':'))

# Create unique model matrix

fireonly_conmod <- lm(formula = sp_rich ~ fire_cat, data = sum_dat,x=T)$x

umm_fireonly <- unique(fireonly_conmod)
rownames(umm_fireonly) <- 1:nrow(umm_fireonly)

# WARNING NUMERIC SUBSETS - put them in the natural order: Intercept (Unburnt), Medium, Burnt

umm_fireonly <- umm_fireonly[c(1,3,2),]
rownames(umm_fireonly) <- 1:nrow(umm_fireonly)

# Create a difference matrix

# Each row must be a vector with a length equal to the number of rows in the unique model matrix (umm_fireonly), e.g. three rows in umm_fireonly matrix will give 3 contrasts. Each row will specify one contrast.

diffm_fireonly <- rbind(
  c(-1,1,0),
  c(-1,0,1),
  c(0,-1,1)
)

# Now we have a unique model matrix
umm_fireonly

# and we have a difference matrix - -1 is baseline, 1 is comparison level to umm
diffm_fireonly

# and we have the names for the contrasts
fireonly_contrast

# fire + location contrasts

# There are 6 contrasts for three categories within each site (hu:hm, hu:hb, hm:hb, pu:pm, pu:pb, pm:pb)

summary(m3_b) # glm.nb model

fireloc_contrast<-data.frame(location=rep(c('Hincks','Pinks'),rep(3,2)),contrast=paste(combn(fireonly.pr$fire_cat,2)[1,],combn(fireonly.pr$fire_cat,2)[2,],sep=':'))

# Create unique model matrix

fireloc_conmod <- lm(abund_25 ~ fire_cat + location, data = sum_dat,x=T)$x

umm_fireloc <- unique(fireloc_conmod)
rownames(umm_fireloc) <- 1:nrow(umm_fireloc)

# WARNING NUMERIC SUBSETS - put them in the natural order: Intercept (Unburnt), Medium, Burnt

umm_fireloc <- umm_fireloc[c(1,3,2,5,6,4),]
rownames(umm_fireloc) <- 1:nrow(umm_fireloc)

# Create a difference matrix

# Each row must be a vector with a length equal to the number of rows in the unique model matrix (umm_fireonly), e.g. three rows in umm_fireonly matrix will give 3 contrasts. Each row will specify one contrast.

diffm_fireloc <- rbind(
  c(-1,1,0,0,0,0),
  c(-1,0,1,0,0,0),
  c(0,-1,1,0,0,0),
  c(0,0,0,-1,1,0),
  c(0,0,0,-1,0,1),
  c(0,0,0,0,-1,1)
)

# Now we have a unique model matrix
umm_fireloc

# and we have a difference matrix - -1 is baseline, 1 is comparison level to umm
diffm_fireloc

# and we have the names for the contrasts
fireloc_contrast

# save.image("04_workspaces/analysed_data.RData")

# ----

# Calculate differences and CI's ----

summary(m1_c) # Species richness, fire only, glm.nb model

# significance based on 'diff', 0 insignificant, 1 significant / lci and uci also show significance based on whether lci + uci cross 0

m1_c_diff<-data.frame(contrast=fireonly_contrast,diff.est(model = m1_c,unique.mod.mat = umm_fireonly,diff.matrix = diffm_fireonly))

m1_c_diff$diff <- ifelse(sign(m1_c_diff$lci)==sign(m1_c_diff$uci),1,0)

summary(m2_c) # Simpsons Index, fire only, Gamma glm

m2_c_diff<-data.frame(contrast=fireonly_contrast,diff.est(model = m2_c,unique.mod.mat = umm_fireonly,diff.matrix = diffm_fireonly))

m2_c_diff$diff <- ifelse(sign(m2_c_diff$lci)==sign(m2_c_diff$uci),1,0)

summary(m3_b) # Abund lowest 25%, fire + location, glm.nb model

m3_b_diff<-data.frame(contrast=fireloc_contrast,diff.est(model = m3_b,unique.mod.mat = umm_fireloc,diff.matrix = diffm_fireloc))

m3_b_diff$diff <- ifelse(sign(m3_b_diff$lci)==sign(m3_b_diff$uci),1,0)

summary (m4_c) # Abund 5% of max, fire only, glm.nb model

m4_c_diff<-data.frame(contrast=fireonly_contrast,diff.est(model = m4_c,unique.mod.mat = umm_fireonly,diff.matrix = diffm_fireonly))

m4_c_diff$diff <- ifelse(sign(m4_c_diff$lci)==sign(m4_c_diff$uci),1,0)

summary (m5_b) # Shannon's index, fire + location, Gamma glm

m5_b_diff<-data.frame(contrast=fireloc_contrast,diff.est(model = m5_b,unique.mod.mat = umm_fireloc,diff.matrix = diffm_fireloc))

m5_b_diff$diff <- ifelse(sign(m5_b_diff$lci)==sign(m5_b_diff$uci),1,0)

summary (m6_b_beta) # Evenness, fire + location, beta gam

m6_b_diff<-data.frame(contrast=fireloc_contrast,diff.est(model = m6_b_beta,unique.mod.mat = umm_fireloc,diff.matrix = diffm_fireloc))

m6_b_diff$diff <- ifelse(sign(m6_b_diff$lci)==sign(m6_b_diff$uci),1,0)

summary (m7_c) # Rare species lowest 25%, fire only, glm.nb

m7_c_diff<-data.frame(contrast=fireonly_contrast,diff.est(model = m7_c,unique.mod.mat = umm_fireonly,diff.matrix = diffm_fireonly))

m7_c_diff$diff <- ifelse(sign(m7_c_diff$lci)==sign(m7_c_diff$uci),1,0)

summary (m8_c) # Rare species 5% max, fire only, glm.nb

m8_c_diff<-data.frame(contrast=fireonly_contrast,diff.est(model = m8_c,unique.mod.mat = umm_fireonly,diff.matrix = diffm_fireonly))

m8_c_diff$diff <- ifelse(sign(m8_c_diff$lci)==sign(m8_c_diff$uci),1,0)

summary (m9_c) # Berger-Parker, fire only, beta gam

m9_c_diff<-data.frame(contrast=fireonly_contrast,diff.est(model = m9_c,unique.mod.mat = umm_fireonly,diff.matrix = diffm_fireonly))

m9_c_diff$diff <- ifelse(sign(m9_c_diff$lci)==sign(m9_c_diff$uci),1,0)

summary (m10_b) # Fisher's alpha, fire + location, gamma

m10_b_diff<-data.frame(contrast=fireloc_contrast,diff.est(model = m10_b,unique.mod.mat = umm_fireloc,diff.matrix = diffm_fireloc))

m10_b_diff$diff <- ifelse(sign(m10_b_diff$lci)==sign(m10_b_diff$uci),1,0)

# ----

# save.image("04_workspaces/analysed_data.RData")

# Summary of all contrasts ----

summary(m1_c) # Species richness, fire only, glm.nb model,
m1_c_diff

summary(m2_c) # Simpsons Index, fire only, Gamma glm
m2_c_diff

summary(m3_b) # Abund lowest 25%, fire + location, glm.nb model
m3_b_diff

summary (m4_c) # Abund 5% of max, fire only, glm.nb model
m4_c_diff

summary (m5_b) # Shannon's index, fire + location, Gamma glm
m5_b_diff

summary (m6_b_beta) # Evenness, fire + location, beta gam
m6_b_diff

summary (m7_c) # Rare species richness, lowest 25%, fire only, glm.nb
m7_c_diff

summary (m8_c) # Rare species richness, 5% max, fire only, glm.nb
m8_c_diff

summary (m9_c) # Berger-Parker, fire only, beta gam
m9_c_diff

summary (m10_b) # Fisher's alpha, fire + location, gamma
m10_b_diff

# ----

# save.image("04_workspaces/analysed_data.RData")

#### PLOTS

#### RARE SPECIES assemblage ----

#### Species richness of rare species on the same figure as abundance of rare species:

# Lowest 25% (14 species). This metric classified species as rare if their abundances fell into the lowest 25% of the assemblage. For these 14 species, we modelled both the richness of this assemblage and the total abundance of all 14 species. 

# 5 % of max abundance
m7_c.pr2
m8_c.pr2
head(sum_dat,3); dim(sum_dat)
# m8_c: sum_dat$sr_5
# m7_c: sum_dat$sr_25
# m4_c: sum_dat$abund_5
# m3_b: sum_dat$abund_25

# Update new axis labels
# make sure they're in the right order (m8_c.pr2$fire_cat;m3_b.pr2$fire_cat)

new.xlab<-c("Long","Medium","Recent")

dev.new(width=7, height=6, dpi=80, pointsize=14, noRStudioGD = T)
par(mfrow=c(2,2), mar=c(4.5,4,1,1), mgp=c(2.2,0.8,0), oma=c(0,0,1,6))

# species richness (5% of max.)
plot(c(1:3),m8_c.pr2$fit, xlim=c(0.5,3.5), pch=20, xaxt="n",ylim=c(0,10),ylab="Richness (5% of max.)",xlab="", las = 1, cex = 2.5, type="n")
arrows(c(1:3),m8_c.pr2$lci85,c(1:3),m8_c.pr2$uci85,length=0,lwd=4,col="grey20")
arrows(c(1:3),m8_c.pr2$lci,c(1:3),m8_c.pr2$uci,length=0.03,code=3,angle=90)
points(c(1:3),m8_c.pr2$fit,cex = 2.5,pch=20)
axis(1,at=c(1:3),labels=F)

axis(1,at=c(1,2,3.15), cex.axis=1,labels=new.xlab,tick=F)
title(mgp=c(2.3,0.8,0),xlab="Fire Category")
m8.tab2
mtext(as.expression(bquote(Delta~"AICc ="~.(paste(round(m8.tab2$Delta_AICc[m8.tab2$Modnames=="fire"]-m1.tab2$Delta_AICc[m8.tab2$Modnames=="null"],2),sep="")))),side=3,line=0.1,adj=1,col="red", cex=0.75)
mtext(text="(a)", side = 3, line = 0.5, adj = 0, cex = 1)
m8_c_diff
# text(x=1:3, y=max(m8_c.pr2$uci)+1,labels=c(letters[1],paste(letters[1],letters[2],sep="",collapse=""),letters[2]))

# plot raw data:
sr5_ub<-sum_dat$sr_5[sum_dat$fire_cat=="Unburnt"]
sr5_m<-sum_dat$sr_5[sum_dat$fire_cat=="Medium"]
sr5_b<-sum_dat$sr_5[sum_dat$fire_cat=="Burnt"]
points(jitter(rep(1,length(sr5_ub)),factor=4),sr5_ub, pch=20,col=rgb(0,0,0,0.2))
points(jitter(rep(2,length(sr5_m)),factor=4),sr5_m, pch=20,col=rgb(0,0,0,0.2))
points(jitter(rep(3,length(sr5_b)),factor=4),sr5_b, pch=20,col=rgb(0,0,0,0.2))

# species richness (Lowest 25%)
plot(c(1:3),m7_c.pr2$fit, xlim=c(0.5,3.5), pch=20, xaxt="n",ylim= c((min(sum_dat$sr_25)),max(m7_c.pr2$uci)+1),ylab="Richness (lowest 25%)",xlab="", las = 1, cex = 2.5, type="n")
arrows(c(1:3),m7_c.pr2$lci85,c(1:3),m7_c.pr2$uci85,length=0,lwd=4,col="grey20")
arrows(c(1:3),m7_c.pr2$lci,c(1:3),m7_c.pr2$uci,length=0.03,code=3,angle=90)
points(c(1:3),m7_c.pr2$fit, pch=20, cex=2.5)
axis(1,at=c(1:3),labels=F)
axis(1,at=c(1,2,3.15),labels=new.xlab,tick=F, cex.axis=1)
title(mgp=c(2.3,0.8,0),xlab="Fire Category")
m7.tab2
mtext(as.expression(bquote(Delta~"AICc ="~.(paste(round(m7.tab2$Delta_AICc[m7.tab2$Modnames=="fire"]-m7.tab2$Delta_AICc[m7.tab2$Modnames=="null"],2),sep="")))), side=3,line=0.1,adj=1,col="red", cex=0.75)
mtext(text="(b)", side = 3, line = 0.5, adj = 0, cex = 1)
m7_c_diff # no differences
# text(x=1:3, y=max(m7_c.pr2$uci)+0.5,labels=rep(letters[1],3))

# plot raw data:
sr25_ub<-sum_dat$sr_25[sum_dat$fire_cat=="Unburnt"]
sr25_m<-sum_dat$sr_25[sum_dat$fire_cat=="Medium"]
sr25_b<-sum_dat$sr_25[sum_dat$fire_cat=="Burnt"]
points(jitter(rep(1,length(sr25_ub)),factor=4),sr25_ub, pch=20,col=rgb(0,0,0,0.2))
points(jitter(rep(2,length(sr25_m)),factor=4),sr25_m, pch=20,col=rgb(0,0,0,0.2))
points(jitter(rep(3,length(sr25_b)),factor=4),sr25_b, pch=20,col=rgb(0,0,0,0.2))

# plot legend:
par(xpd=NA)
legend(x=4,y=max(m7_c.pr2$uci)+1.1, title = "Location", legend = c("Fire only", "Hincks","Pinkaw."), pt.cex = 1.5, pch = c(16, 15, 17), bty = "n", title.adj=0)
par(xpd=F)

# Abundance (5% of max.)
plot(c(1:3),m4_c.pr2$fit, xlim=c(0.5,3.5), pch=20, xaxt="n",ylim= c(0,max(m4_c.pr2$uci)+2),ylab="Abundance (5% of max.)",xlab="", las = 1, cex = 2.5, type="n")
arrows(c(1:3),m4_c.pr2$lci85,c(1:3),m4_c.pr2$uci85,length=0,lwd=4,col="grey20")
arrows(c(1:3),m4_c.pr2$lci,c(1:3),m4_c.pr2$uci,length=0.03,code=3,angle=90)
points(c(1:3),m4_c.pr2$fit, cex=2.5, pch=20)
axis(1,at=c(1:3),labels=F)
axis(1,at=c(1,2,3.15),cex.axis=1,labels=new.xlab,tick=F)
title(mgp=c(2.3,0.8,0),xlab="Fire Category")
mtext(as.expression(bquote(Delta~"AICc ="~.(paste(round(m4.tab2$Delta_AICc[m4.tab2$Modnames=="fire"]-m4.tab2$Delta_AICc[m4.tab2$Modnames=="null"],2),sep="")))), side=3,line=0.1,adj=1,col="red", cex=0.75)
mtext(text="(c)", side = 3, line = 0.5, adj = 0, cex = 1)
m4_c_diff # no differences
# text(x=1:3, y=max(m4_c.pr2$uci)+1,labels=rep(letters[1],3))

# plot raw data:
abund5_ub<-sum_dat$abund_5[sum_dat$fire_cat=="Unburnt"]
abund5_m<-sum_dat$abund_5[sum_dat$fire_cat=="Medium"]
abund5_b<-sum_dat$abund_5[sum_dat$fire_cat=="Burnt"]

points(jitter(rep(1,length(abund5_ub)),factor=4),abund5_ub, pch=20,col=rgb(0,0,0,0.2))
points(jitter(rep(2,length(abund5_m)),factor=4),abund5_m, pch=20,col=rgb(0,0,0,0.2))
points(jitter(rep(3,length(abund5_b)),factor=4),abund5_b, pch=20,col=rgb(0,0,0,0.2))

# Abundance (Lowest 25%)
plot(c(1:3)-0.15,m3_b.pr2$fit[m3_b.pr2$location=="Hincks"], xlim=c(0.5,3.5), pch=15, xaxt="n",ylim= c((min(sum_dat$abund_25)),max(m3_b.pr2$uci)),ylab="Abundance (lowest 25%)",xlab="", las = 1, cex = 1.5, type="n")
arrows(c(1:3)-0.15,m3_b.pr2$lci85[m3_b.pr2$location=="Hincks"],c(1:3)-0.15,m3_b.pr2$uci85[m3_b.pr2$location=="Hincks"],length=0,lwd=4,col="grey20")
arrows(c(1:3)+0.15,m3_b.pr2$lci85[m3_b.pr2$location=="Pinks"],c(1:3)+0.15,m3_b.pr2$uci85[m3_b.pr2$location=="Pinks"],length=0,lwd=4,col="grey20")
arrows(c(1:3)-0.15,m3_b.pr2$lci[m3_b.pr2$location=="Hincks"],c(1:3)-0.15,m3_b.pr2$uci[m3_b.pr2$location=="Hincks"],length=0.03,code=3,angle=90)
arrows(c(1:3)+0.15,m3_b.pr2$lci[m3_b.pr2$location=="Pinks"],c(1:3)+0.15,m3_b.pr2$uci[m3_b.pr2$location=="Pinks"],length=0.03,code=3,angle=90)
points(c(1:3)-0.15,m3_b.pr2$fit[m3_b.pr2$location=="Hincks"], xlim=c(0.5,3.5), pch=15, cex = 1.5)
points(c(1:3)+0.15,m3_b.pr2$fit[m3_b.pr2$location=="Pinks"], xlim=c(0.5,3.5), pch=17, cex = 1.5)

axis(1,at=c(1:3),labels=F)
axis(1,at=c(1,2,3.15), cex.axis=1,labels=new.xlab,tick=F)
title(mgp=c(2.3,0.8,0),xlab="Fire Category")
mtext(as.expression(bquote(Delta~"AICc ="~.(paste(round(m3.tab2$Delta_AICc[m3.tab2$Modnames=="location"]-m3.tab2$Delta_AICc[m3.tab2$Modnames=="null"],2),sep="")))), side=3,line=0.1,adj=1,col="red", cex=0.75)
mtext(text="(d)", side = 3, line = 0.5, adj = 0, cex = 1)
m3_b_diff
# text(x=1:3, y=max(m3_b.pr2$uci)+1,labels=c(letters[1],paste(letters[1],letters[2],sep="",collapse=""),letters[2]))

# plot raw data:
abund25_ubH<-sum_dat$abund_25[sum_dat$fire_cat=="Unburnt" & sum_dat$location=="Hincks"]
abund25_ubP<-sum_dat$abund_25[sum_dat$fire_cat=="Unburnt" & sum_dat$location=="Pinks"]
abund25_mH<-sum_dat$abund_25[sum_dat$fire_cat=="Medium" & sum_dat$location=="Hincks"]
abund25_mP<-sum_dat$abund_25[sum_dat$fire_cat=="Medium" & sum_dat$location=="Pinks"]
abund25_bH<-sum_dat$abund_25[sum_dat$fire_cat=="Burnt" & sum_dat$location=="Hincks"]
abund25_bP<-sum_dat$abund_25[sum_dat$fire_cat=="Burnt" & sum_dat$location=="Pinks"]

points(jitter(rep(1-0.15,length(abund25_ubH)),factor=4),abund25_ubH, pch=15,col=rgb(0,0,0,0.2))
points(jitter(rep(1+0.15,length(abund25_ubP)),factor=4),abund25_ubP, pch=17,col=rgb(0,0,0,0.2))

points(jitter(rep(2-0.15,length(abund25_mH)),factor=4),abund25_mH, pch=15,col=rgb(0,0,0,0.2))
points(jitter(rep(2+0.15,length(abund25_mP)),factor=4),abund25_mP, pch=17,col=rgb(0,0,0,0.2))

points(jitter(rep(3-0.15,length(abund25_bH)),factor=4),abund25_bH, pch=15,col=rgb(0,0,0,0.2))
points(jitter(rep(3+0.15,length(abund25_bP)),factor=4),abund25_bP, pch=17,col=rgb(0,0,0,0.2))

# ----

#### ALL SPECIES assemblage ----

#### Species richness all species:

# Remove the a, b, c, etc. from the plots 
# update the whole community plot to include the bp and fa
# Presented in the MS as follows: species richness, Shannon-Weiner’s Diversity Index (hereafter Shannon’s Diversity), Inverse Simpson’s Diversity Index (hereafter Simpson’s Diversity), Fisher’s alpha, Berger-Parker Index, and evenness

# m1_c: sum_dat$sp_rich (fire only); m1.tab
# m5_b: sum_dat$shann_ind (fire + location); m5.tab
# m2_c: sum_dat$simps_ind2 (fire only); m2.tab
# m10_b: sum_dat$fa.all (fire + location); m10.tab
# m9_c: sum_dat$bp (fire only); m9.tab
# m6_b_beta: sum_dat$even2 (fire + location); m6.tab

dev.new(width=7, height=9, dpi=80, pointsize=17.5, noRStudioGD = T)
par(mfrow=c(3,2), mar=c(4.5,4,1,1), mgp=c(2.3,0.6,0), oma=c(0,0,1,6))

# species richness; fire only

m1.ciraw<-c(min(sum_dat$sp_rich),max(sum_dat$sp_rich),min(m1_c.pr2$lci),max(m1_c.pr2$uci))

plot(c(1:3),m1_c.pr2$fit, xlim=c(0.5,3.5), pch=20, xaxt="n",ylim=c(min(m1.ciraw),max(m1.ciraw)),ylab="Species Richness",xlab="", las = 1, cex = 2.5, type="n")
arrows(c(1:3),m1_c.pr2$lci85,c(1:3),m1_c.pr2$uci85,length=0,lwd=4,col="grey20")
arrows(c(1:3),m1_c.pr2$lci,c(1:3),m1_c.pr2$uci,length=0.03,code=3,angle=90)
points(c(1:3),m1_c.pr2$fit,cex = 2.5,pch=20)
axis(1,at=c(1:3),labels=F)
axis(1,at=c(1,2,3.15), cex.axis=1,labels=new.xlab,tick=F)
title(mgp=c(2.1,0.8,0),xlab="Fire Category")
m1.tab2
mtext(as.expression(bquote(Delta~"AICc ="~.(paste(round(m1.tab2$Delta_AICc[m1.tab2$Modnames=="fire"]-m1.tab2$Delta_AICc[m1.tab2$Modnames=="null"],2),sep="")))), side=3,line=0.1,adj=1,col="red", cex=0.70)
mtext(text="(a)", side = 3, line = 0.5, adj = 0, cex = 0.8)
m1_c_diff # no differences
# text(x=1:3, y=max(m1_c.pr2$uci)+2,labels=rep(letters[1],3))

# plot raw data:
sr_ub<-sum_dat$sp_rich[sum_dat$fire_cat=="Unburnt"]
sr_m<-sum_dat$sp_rich[sum_dat$fire_cat=="Medium"]
sr_b<-sum_dat$sp_rich[sum_dat$fire_cat=="Burnt"]

points(jitter(rep(1,length(sr_ub)),factor=4),sr_ub, pch=20,col=rgb(0,0,0,0.2))
points(jitter(rep(2,length(sr_m)),factor=4),sr_m, pch=20,col=rgb(0,0,0,0.2))
points(jitter(rep(3,length(sr_b)),factor=4),sr_b, pch=20,col=rgb(0,0,0,0.2))

# shann_ind; fire+location

plot(c(1:3)-0.15,m5_b.pr2$fit[m5_b.pr2$location=="Hincks"], xlim=c(0.5,3.5), pch=15, xaxt="n",ylim= c((min(m5_b.pr2$lci)),max(m5_b.pr2$uci)),ylab="Shannon's Diversity",xlab="", las = 1, cex = 1.5, type="n")
arrows(c(1:3)-0.15,m5_b.pr2$lci85[m5_b.pr2$location=="Hincks"],c(1:3)-0.15,m5_b.pr2$uci85[m5_b.pr2$location=="Hincks"],length=0,lwd=4,col="grey20")
arrows(c(1:3)+0.15,m5_b.pr2$lci85[m5_b.pr2$location=="Pinks"],c(1:3)+0.15,m5_b.pr2$uci85[m5_b.pr2$location=="Pinks"],length=0,lwd=4,col="grey20")
points(c(1:3)-0.15,m5_b.pr2$fit[m5_b.pr2$location=="Hincks"], xlim=c(0.5,3.5), pch=15, cex = 1.5)
points(c(1:3)+0.15,m5_b.pr2$fit[m5_b.pr2$location=="Pinks"], xlim=c(0.5,3.5), pch=17, cex = 1.5)
arrows(c(1:3)-0.15,m5_b.pr2$lci[m5_b.pr2$location=="Hincks"],c(1:3)-0.15,m5_b.pr2$uci[m5_b.pr2$location=="Hincks"],length=0.03,code=3,angle=90)
arrows(c(1:3)+0.15,m5_b.pr2$lci[m5_b.pr2$location=="Pinks"],c(1:3)+0.15,m5_b.pr2$uci[m5_b.pr2$location=="Pinks"],length=0.03,code=3,angle=90)
axis(1,at=c(1:3),labels=F)
axis(1,at=c(1,2,3.15),labels=new.xlab,tick=F)
title(mgp=c(2.1,0.8,0),xlab="Fire Category")
mtext(as.expression(bquote(Delta~"AICc ="~.(paste(round(m5.tab2$Delta_AICc[m5.tab2$Modnames=="location"]-m5.tab2$Delta_AICc[m5.tab2$Modnames=="null"],2),sep="")))), side=3,line=0.1,adj=1,col="darkgreen", cex=0.7)
mtext(text="(b)", side = 3, line = 0.5, adj = 0, cex = 0.8)
m5_b_diff
# text(x=1:3, y=max(m5_b.pr2$uci)+0.2,labels=c(letters[1],rep(letters[2],2)))

# plot raw data:
sh_ubH<-sum_dat$shann_ind[sum_dat$fire_cat=="Unburnt" & sum_dat$location=="Hincks"]
sh_ubP<-sum_dat$shann_ind[sum_dat$fire_cat=="Unburnt" & sum_dat$location=="Pinks"]
sh_mH<-sum_dat$shann_ind[sum_dat$fire_cat=="Medium" & sum_dat$location=="Hincks"]
sh_mP<-sum_dat$shann_ind[sum_dat$fire_cat=="Medium" & sum_dat$location=="Pinks"]
sh_bH<-sum_dat$shann_ind[sum_dat$fire_cat=="Burnt" & sum_dat$location=="Hincks"]
sh_bP<-sum_dat$shann_ind[sum_dat$fire_cat=="Burnt" & sum_dat$location=="Pinks"]

points(jitter(rep(1-0.15,length(sh_ubH)),factor=4),sh_ubH, pch=15,col=rgb(0,0,0,0.2))
points(jitter(rep(1+0.15,length(sh_ubP)),factor=4),sh_ubP, pch=17,col=rgb(0,0,0,0.2))

points(jitter(rep(2-0.15,length(sh_mH)),factor=4),sh_mH, pch=15,col=rgb(0,0,0,0.2))
points(jitter(rep(2+0.15,length(sh_mP)),factor=4),sh_mP, pch=17,col=rgb(0,0,0,0.2))

points(jitter(rep(3-0.15,length(sh_bH)),factor=4),sh_bH, pch=15,col=rgb(0,0,0,0.2))
points(jitter(rep(3+0.15,length(sh_bP)),factor=4),sh_bP, pch=17,col=rgb(0,0,0,0.2))

# plot legend:
par(xpd=NA)
legend(x=4,y=max(m5_b.pr2$uci), title = "Location", legend = c("Fire only", "Hincks","Pinkaw."), pt.cex = 1.5, pch = c(16, 15, 17), bty = "n", title.adj=0)
par(xpd=F)

# simps diversity index; fire only

# plot raw data:
# run first as this is used for ylim
si2_ub<-sum_dat$simps_ind2[sum_dat$fire_cat=="Unburnt"]
si2_m<-sum_dat$simps_ind2[sum_dat$fire_cat=="Medium"]
si2_b<-sum_dat$simps_ind2[sum_dat$fire_cat=="Burnt"]

plot(c(1:3),m2_c.pr2$fit, xlim=c(0.5,3.5), pch=20, xaxt="n",ylim=c(min(si2_b),max(sum_dat$simps_ind2)+1),ylab="Simpson's Diversity",xlab="", las = 1, cex = 2.5, type="n")
arrows(c(1:3),m2_c.pr2$lci85,c(1:3),m2_c.pr2$uci85,length=0,lwd=4,col="grey20")
arrows(c(1:3),m2_c.pr2$lci,c(1:3),m2_c.pr2$uci,length=0.03,code=3,angle=90)
points(c(1:3),m2_c.pr2$fit, cex=2.5, pch=20)
axis(1,at=c(1:3),labels=F)
axis(1,at=c(1,2,3.15), cex.axis=1,labels=new.xlab,tick=F)
title(mgp=c(2.1,0.8,0),xlab="Fire Category")
mtext(as.expression(bquote(Delta~"AICc ="~.(paste(round(m2.tab2$Delta_AICc[m2.tab2$Modnames=="fire"]-m2.tab2$Delta_AICc[m2.tab2$Modnames=="null"],2),sep="")))), side=3,line=0.1,adj=1,col="red", cex=0.7)
mtext(text="(c)", side = 3, line = 0.5, adj = 0, cex = 0.8)
m2_c_diff
# text(x=1:3, y=max(m2_c.pr2$uci)+1.5,labels=c(letters[1],rep(letters[2],2)))

points(jitter(rep(1,length(si2_ub)),factor=4),si2_ub, pch=20,col=rgb(0,0,0,0.2))
points(jitter(rep(2,length(si2_m)),factor=4),si2_m, pch=20,col=rgb(0,0,0,0.2))
points(jitter(rep(3,length(si2_b)),factor=4),si2_b, pch=20,col=rgb(0,0,0,0.2))

# Fisher's alpha; fire+location

# m10_b: sum_dat$fa.all (fire + location); m10.tab
# m10_b.pr2

plot(c(1:3)-0.15,m10_b.pr2$fit[m10_b.pr2$location=="Hincks"], xlim=c(0.5,3.5), pch=15, xaxt="n",ylim= c((min(sum_dat$fa.all)),max(m10_b.pr2$uci)),ylab="",xlab="", las = 1, cex = 1.5, type="n")

title(ylab=as.expression(bquote("Fisher's"~alpha)),mgp=c(2.2,0.8,0))

arrows(c(1:3)-0.15,m10_b.pr2$lci85[m10_b.pr2$location=="Hincks"],c(1:3)-0.15,m10_b.pr2$uci85[m10_b.pr2$location=="Hincks"],length=0,lwd=4,col="grey20")
arrows(c(1:3)+0.15,m10_b.pr2$lci85[m10_b.pr2$location=="Pinks"],c(1:3)+0.15,m10_b.pr2$uci85[m10_b.pr2$location=="Pinks"],length=0,lwd=4,col="grey20")

arrows(c(1:3)-0.15,m10_b.pr2$lci[m10_b.pr2$location=="Hincks"],c(1:3)-0.15,m10_b.pr2$uci[m10_b.pr2$location=="Hincks"],length=0.03,code=3,angle=90)
arrows(c(1:3)+0.15,m10_b.pr2$lci[m10_b.pr2$location=="Pinks"],c(1:3)+0.15,m10_b.pr2$uci[m10_b.pr2$location=="Pinks"],length=0.03,code=3,angle=90)

points(c(1:3)-0.15,m10_b.pr2$fit[m10_b.pr2$location=="Hincks"], xlim=c(0.5,3.5), pch=15, cex = 1.5)
points(c(1:3)+0.15,m10_b.pr2$fit[m10_b.pr2$location=="Pinks"], xlim=c(0.5,3.5), pch=17, cex = 1.5)

axis(1,at=c(1:3),labels=F)
axis(1,at=c(1,2,3.15),labels=new.xlab,tick=F)
title(mgp=c(2.1,0.8,0),xlab="Fire Category")
mtext(as.expression(bquote(Delta~"AICc ="~.(paste(round(m10.tab$Delta_AICc[m10.tab$Modnames=="location"]-m10.tab$Delta_AICc[m10.tab$Modnames=="null"],2),sep="")))), side=3,line=0.1,adj=1,col="darkgreen", cex=0.7)
mtext(text="(d)", side = 3, line = 0.5, adj = 0, cex = 0.8)

# plot raw data:
fa_ubH<-sum_dat$fa.all[sum_dat$fire_cat=="Unburnt" & sum_dat$location=="Hincks"]
fa_ubP<-sum_dat$fa.all[sum_dat$fire_cat=="Unburnt" & sum_dat$location=="Pinks"]
fa_mH<-sum_dat$fa.all[sum_dat$fire_cat=="Medium" & sum_dat$location=="Hincks"]
fa_mP<-sum_dat$fa.all[sum_dat$fire_cat=="Medium" & sum_dat$location=="Pinks"]
fa_bH<-sum_dat$fa.all[sum_dat$fire_cat=="Burnt" & sum_dat$location=="Hincks"]
fa_bP<-sum_dat$fa.all[sum_dat$fire_cat=="Burnt" & sum_dat$location=="Pinks"]

points(jitter(rep(1-0.15,length(fa_ubH)),factor=4),fa_ubH, pch=15,col=rgb(0,0,0,0.2))
points(jitter(rep(1+0.15,length(fa_ubP)),factor=4),fa_ubP, pch=17,col=rgb(0,0,0,0.2))

points(jitter(rep(2-0.15,length(fa_mH)),factor=4),fa_mH, pch=15,col=rgb(0,0,0,0.2))
points(jitter(rep(2+0.15,length(fa_mP)),factor=4),fa_mP, pch=17,col=rgb(0,0,0,0.2))

points(jitter(rep(3-0.15,length(fa_bH)),factor=4),fa_bH, pch=15,col=rgb(0,0,0,0.2))
points(jitter(rep(3+0.15,length(fa_bP)),factor=4),fa_bP, pch=17,col=rgb(0,0,0,0.2))

# Berger-Parker; fire only

# m9_c: sum_dat$bp (fire only); m9.tab
# m9_c.pr2

m9.ciraw<-c(min(sum_dat$bp_ind),max(sum_dat$bp_ind),min(m9_c.pr2$lci),max(m9_c.pr2$uci))

plot(c(1:3),m9_c.pr2$fit, xlim=c(0.5,3.5), pch=20, xaxt="n",ylim=c(min(m9.ciraw),max(m9.ciraw)),ylab="Berger Parker Index",xlab="", las = 1, cex = 2.5, type="n")

arrows(c(1:3),m9_c.pr2$lci85,c(1:3),m9_c.pr2$uci85,length=0,lwd=4,col="grey20")

arrows(c(1:3),m9_c.pr2$lci,c(1:3),m9_c.pr2$uci,length=0.03,code=3,angle=90)
points(c(1:3),m9_c.pr2$fit,cex=2.5,pch=20)

axis(1,at=c(1:3),labels=F)
axis(1,at=c(1,2,3.15), cex.axis=1,labels=new.xlab,tick=F)
title(mgp=c(2.1,0.8,0),xlab="Fire Category")
mtext(as.expression(bquote(Delta~"AICc ="~.(paste(round(m9.tab$Delta_AICc[m9.tab$Modnames=="fire"]-m9.tab$Delta_AICc[m9.tab$Modnames=="null"],2),sep="")))), side=3,line=0.1,adj=1,col="red", cex=0.7)
mtext(text="(e)", side = 3, line = 0.5, adj = 0, cex = 0.8)

# plot raw data:
bp_ub<-sum_dat$bp_ind[sum_dat$fire_cat=="Unburnt"]
bp_m<-sum_dat$bp_ind[sum_dat$fire_cat=="Medium"]
bp_b<-sum_dat$bp_ind[sum_dat$fire_cat=="Burnt"]

points(jitter(rep(1,length(bp_ub)),factor=4),bp_ub, pch=20,col=rgb(0,0,0,0.2))
points(jitter(rep(2,length(bp_m)),factor=4),bp_m, pch=20,col=rgb(0,0,0,0.2))
points(jitter(rep(3,length(bp_b)),factor=4),bp_b, pch=20,col=rgb(0,0,0,0.2))

# evenness

# m6_b_beta: sum_dat$even2 (fire + location); m6.tab
# m6_b.pr2

m6.ciraw<-c(min(sum_dat$even2),max(sum_dat$even2),min(m6_b.pr2$lci),max(m6_b.pr2$uci))

plot(c(1:3)-0.15,m6_b.pr2$fit[m6_b.pr2$location=="Hincks"], xlim=c(0.5,3.5), pch=15, xaxt="n",ylim=c((min(m6.ciraw)),max(m6.ciraw)),ylab="Evenness",xlab="", las = 1, cex = 1.5, type="n")

arrows(c(1:3)-0.15,m6_b.pr2$lci85[m6_b.pr2$location=="Hincks"],c(1:3)-0.15,m6_b.pr2$uci85[m6_b.pr2$location=="Hincks"],length=0,lwd=4,col="grey20")
arrows(c(1:3)+0.15,m6_b.pr2$lci85[m6_b.pr2$location=="Pinks"],c(1:3)+0.15,m6_b.pr2$uci85[m6_b.pr2$location=="Pinks"],length=0,lwd=4,col="grey20")

arrows(c(1:3)-0.15,m6_b.pr2$lci[m6_b.pr2$location=="Hincks"],c(1:3)-0.15,m6_b.pr2$uci[m6_b.pr2$location=="Hincks"],length=0.03,code=3,angle=90)
arrows(c(1:3)+0.15,m6_b.pr2$lci[m6_b.pr2$location=="Pinks"],c(1:3)+0.15,m6_b.pr2$uci[m6_b.pr2$location=="Pinks"],length=0.03,code=3,angle=90)

points(c(1:3)-0.15,m6_b.pr2$fit[m6_b.pr2$location=="Hincks"], xlim=c(0.5,3.5), pch=15, cex = 1.5)
points(c(1:3)+0.15,m6_b.pr2$fit[m6_b.pr2$location=="Pinks"], xlim=c(0.5,3.5), pch=17, cex = 1.5)

axis(1,at=c(1:3),labels=F)
axis(1,at=c(1,2,3.15),labels=new.xlab,tick=F)
title(mgp=c(2.1,0.8,0),xlab="Fire Category")
mtext(as.expression(bquote(Delta~"AICc ="~.(paste(round(m6.tab$Delta_AICc[m6.tab$Modnames=="location"]-m6.tab$Delta_AICc[m6.tab$Modnames=="null"],2),sep="")))), side=3,line=0.1,adj=1,col="red", cex=0.7)

mtext(text="(f)", side = 3, line = 0.5, adj = 0, cex = 0.8)

# plot raw data:
ev_ubH<-sum_dat$even2[sum_dat$fire_cat=="Unburnt" & sum_dat$location=="Hincks"]
ev_ubP<-sum_dat$even2[sum_dat$fire_cat=="Unburnt" & sum_dat$location=="Pinks"]
ev_mH<-sum_dat$even2[sum_dat$fire_cat=="Medium" & sum_dat$location=="Hincks"]
ev_mP<-sum_dat$even2[sum_dat$fire_cat=="Medium" & sum_dat$location=="Pinks"]
ev_bH<-sum_dat$even2[sum_dat$fire_cat=="Burnt" & sum_dat$location=="Hincks"]
ev_bP<-sum_dat$even2[sum_dat$fire_cat=="Burnt" & sum_dat$location=="Pinks"]

points(jitter(rep(1-0.15,length(ev_ubH)),factor=4),ev_ubH, pch=15,col=rgb(0,0,0,0.2))
points(jitter(rep(1+0.15,length(ev_ubP)),factor=4),ev_ubP, pch=17,col=rgb(0,0,0,0.2))

points(jitter(rep(2-0.15,length(ev_mH)),factor=4),ev_mH, pch=15,col=rgb(0,0,0,0.2))
points(jitter(rep(2+0.15,length(ev_mP)),factor=4),ev_mP, pch=17,col=rgb(0,0,0,0.2))

points(jitter(rep(3-0.15,length(ev_bH)),factor=4),ev_bH, pch=15,col=rgb(0,0,0,0.2))
points(jitter(rep(3+0.15,length(ev_bP)),factor=4),ev_bP, pch=17,col=rgb(0,0,0,0.2))

# ----

# save.image("04_workspaces/analysed_data.RData")

# Rank abundance plot for appendix ----

# Use captures / 1000 trap nights, 2 seasons of data combined
head(sp_div2, 3); dim(sp_div2)

sp_abund<-data.frame(species=names(colSums(sp_div2)), total_abund=colSums(sp_div2))
sp_abund<-sp_abund[order(sp_abund$total_abund),]
rownames(sp_abund)<-1:nrow(sp_abund)
head(sp_abund); dim(sp_abund)
tail(sp_abund)

dev.new(height = 8, width = 12, noRStudioGD = T, dpi=70, pointsize=20)
par(mfrow = c(2,2), mar=c(4,5,2,1), mgp=c(2.3,0.7,0))

hist(sp_abund$total_abund, breaks=seq(0,500,by=10), main = "", ylab="Number of species", xlab="",las=1, col="grey80")
title(xlab="Abundance", mgp=c(2,1,0))
title(main="(a) All Data", mgp=c(2,1,0), adj = 0, font.main = 1, cex.main=0.95)

lab1<-tail(sp_abund,4)
lab1names<-c("Ctenotus euclae","Liopholis inornata","Nephrurus stellatus","Ctenophorus fordi")
text(c(lab1$total_abund[1]-10,lab1$total_abund[2]+10,lab1$total_abund[c(3,4)]),2,font=3,labels=lab1names,srt=90, adj=0, col="darkorange4", cex=0.8)

hist(sp_abund$total_abund[sp_abund$total_abund<100], breaks=seq(0,100,by=10), main="", ylab="Number of species", xlab="",las=1, col="grey80")
title(xlab="Abundance", mgp=c(2,1,0))
title(main="(b) Less than 100", mgp=c(2,1,0), adj = 0, font.main = 1, cex.main=0.95)

lab2<-tail(sp_abund[sp_abund$total_abund<100,],3)
lab2names<-c("Ctenotus atlas","Lerista edwardsae","Ctenophorus cristatus")
text(c(lab2$total_abund[1]+4,lab2$total_abund[c(2,3)]),2,font=3,labels=lab2names,srt=90, adj=0, col="darkorange4", cex=0.8)

hist(sp_abund$total_abund[sp_abund$total_abund<20], breaks=seq(0,20,by=1), main="", ylab="Number of species", xlab="",las=1, col="grey80", adj = 0)
title(xlab="Abundance", mgp=c(2,1,0))
title(main="(c) Less than 20", mgp=c(2,1,0), adj = 0,font.main = 1, cex.main=0.95)

lab3<-tail(sp_abund[sp_abund$total_abund<20,],2)
lab3names<-c("Anilios bituberculatus","Anilios australis")
text(lab3$total_abund-0.5,1.5,font=3,labels=lab3names,srt=90, adj=0, col="darkorange4", cex=0.8)

hist(sp_abund$total_abund[sp_abund$total_abund<10], breaks=seq(0,10,by=1), main="", ylab="Number of species", xlab="",las=1, col="grey80", adj = 0)
title(xlab="Abundance", mgp=c(2,1,0))
title(main="(d) Less than 10", mgp=c(2,1,0), adj = 0, font.main = 1, cex.main=0.95)

# ----

# Coefficients for top models ----

# Whole reptile community:

# m1_c: sum_dat$sp_rich (fire only); m1.tab
# m5_b: sum_dat$shann_ind (fire + location); m5.tab
# m2_c: sum_dat$simps_ind2 (fire only); m2.tab
# m10_b: sum_dat$fa.all (fire + location); m10.tab
# m9_c: sum_dat$bp (fire only); m9.tab
# m6_b_beta: sum_dat$even2 (fire + location); m6.tab

all.coef<-rbind(
  matrix(data=c("Species richness","","",""),nrow=1,ncol=4),
  round(summary(m1_c)$coefficients,3)[,1:4], # Species richness, fire only, glm.nb model
  
  matrix(data=c("Shannon's Diversity Index", "","",""),nrow=1,ncol=4),
  round(summary(m5_b)$coefficients,3)[,1:4], # Shannon's index, fire + location, Gamma glm
  matrix(data=c("Simpson's Diversity Index","","",""),nrow=1,ncol=4),
  round(summary(m2_c)$coefficients,3)[,1:4], # Simpsons Index, fire only, Gamma glm
  
  matrix(data=c("Fisher's Alpha","","",""),nrow=1,ncol=4),
  round(summary(m10_b)$coefficients,3)[,1:4], # Fisher's Alpha, fire + location, Gamma glm
  
  matrix(data=c("Berger-Parker","","",""),nrow=1,ncol=4),
  round(summary(m9_c)$p.table,3)[,1:4], # Berger-Parker, fire only, beta gam
  
  matrix(data=c("Evenness", "","",""),nrow=1,ncol=4),
  round(summary(m6_b_beta)$p.table,3)[,1:4], # Evenness, fire + location, beta gam

# Rare species:

# m8_c: sum_dat$sr_5 (fire only); m8_c.pr2; m8.tab
# m7_c: sum_dat$sr_25 (fire only); m7_c.pr2; m7.tab
# m4_c: sum_dat$abund_5 (fire only); m4_c.pr2; m4.tab
# m3_b: sum_dat$abund_25 (fire + location); m3_b.pr2; m3.tab

matrix(data=c("Richness (5% of max.)", "","",""),nrow=1,ncol=4),
round(summary(m8_c)$coefficients,3)[,1:4], # Rare species richness 5% max, fire only, glm.nb

matrix(data=c("Richness (lowest 25%)", "","",""),nrow=1,ncol=4),
round(summary(m7_c)$coefficients,3)[,1:4], # Rare species richness, lowest 25%, fire only, glm.nb

matrix(data=c("Abundance (5% of max.)", "","",""),nrow=1,ncol=4),
round(summary(m4_c)$coefficients,3)[,1:4], # Abund 5% of max, fire only, glm.nb 

matrix(data=c("Abundance (lowest 25%)", "","",""),nrow=1,ncol=4),
round(summary(m3_b)$coefficients,3)[,1:4] # Abund lowest 25%, fire + location, glm.nb 
)

all.coef2<-as.data.frame(all.coef)
# write.table(all.coef2,file="all.coef2.txt",quote=F,sep="\t",row.names = T)

# ----

# Plot effect sizes for all contrasts (not used in new MS) ---- 

diff.list<-list(
  # Species richness, fire only, glm.nb model
  m1_c_diff, 
  # Shannon's index, fire + location, Gamma glm
  m5_b_diff, 
  # Simpsons Index, fire only, Gamma glm
  m2_c_diff, 
  # Evenness, fire only, Gamma glm
  m6_c_diff, 
  # Rare species richness 5% max, fire only, glm.nb
  m8_c_diff, 
  # Rare species richness, lowest 25%, fire only, glm.nb
  m7_c_diff, 
  # Abund 5% of max, fire only, glm.nb model
  m4_c_diff, 
  # Abund lowest 25%, fire + location, glm.nb model
  m3_b_diff)

diff.lab<-c("Species richness","Shannon's Diversity",  "Simpson's Diversity ","Evenness", "Richness (5% of max.)", "Richness (lowest 25%)", "Abundance (5% of max.)", "Abundance (lowest 25%)")

dev.new(width=4, height=7, dpi=80, pointsize=12, noRStudioGD = T)
par(mfrow=c(4,2), mar=c(4,4,2,1), mgp=c(2,0.8,0))

for(i in 1:length(diff.list)){
  
  list.thisrun<-diff.list[[i]]
  
  if(nrow(list.thisrun)==3){
    plot(rev(list.thisrun$est), 1:3, xlim=c(min(list.thisrun$lci),max(list.thisrun$uci)), yaxt="n", pch=16, ylim=c(0.5,3.5), xlab="Difference ± 95 % CI", ylab="", main="", font.main=1)
    mtext(paste("(",letters[i],") ",diff.lab[i],sep=""),side=3, line=0.5, cex=0.65, adj=0)
    arrows(rev(list.thisrun$lci),1:3,rev(list.thisrun$uci),1:3,length=0, lty=1)
    arrows(0,0,0,4,length=0, lty=2)
    axis(side=2,at=1:3,labels=rev(c("U:M","U:B","M:B")), las=1)
  } # close if nrow == 3
  
  if(nrow(list.thisrun)==6){
    plot(rev(list.thisrun$est), 1:6, xlim=c(min(list.thisrun$lci),max(list.thisrun$uci)), yaxt="n", pch=16, ylim=c(0.5,6.5), xlab="Difference ± 95 % CI", ylab="", main="", font.main=1)
    mtext(paste("(",letters[i],") ",diff.lab[i],sep=""),side=3, line=0.5, cex=0.65, adj=0)
    arrows(rev(list.thisrun$lci),1:6,rev(list.thisrun$uci),1:6,length=0, lty=1)
    arrows(0,0,0,7,length=0, lty=2)
    axis(side=2,at=1:6,labels=rev(c("H U:M","H U:B","H M:B","P U:M","P U:B","P M:B")), las=1)
  } # close if nrow == 6
}

# ----

