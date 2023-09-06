

# 02_analyse_data

#use this if starting fresh
#load("04_workspaces/processed_data.RData")
#already analysed data
load("04_workspaces/analysed_data.RData")
#save.image("04_workspaces/analysed_data.RData")


library("MASS")
library(AICcmodavg)

# Load functions:
invisible(lapply(paste("02_functions/",dir("02_functions"),sep=""), function(x) source(x)))

#data fully processed and ready to process 9th May 2023 in sum_dat

head(sum_dat, 6);dim(sum_dat)

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

m1.set<-list("fire x location"= m1_a, "location"= m1_b, "fire"= m1_c, "null"= m1_d)
m1.tab<-aictab(cand.set = m1.set, second.ord = T, sort = T)

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

m2.set<-list("fire x location"= m2_a, "location"= m2_b, "fire"= m2_c, "null"= m2_d)
m2.tab<-aictab(cand.set = m2.set, second.ord = T, sort = T)

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

m3.set<-list("fire x location"= m3_a, "location"= m3_b, "fire"= m3_c, "null"= m3_d)
m3.tab<-aictab(cand.set = m3.set, second.ord = T, sort = T)

#write.table(m3.tab,file="m3_tab.txt", quote = F, sep = "\t", row.names=F)

#predictions for models to plot graphs
fireonly.pr<-data.frame(fire_cat = factor(levels(sum_dat$fire_cat),levels = levels(sum_dat$fire_cat)))
m3_c.pr<-predict(object = m3_c, newdata = fireonly.pr,type = "response", se.fit = T)
m3_c.pr2<-data.frame(fireonly.pr)
m3_c.pr2$fit<-m3_c.pr$fit
m3_c.pr2$se<-m3_c.pr$se
m3_c.pr2$lci<-m3_c.pr$fit-(m3_c.pr2$se*1.96)
m3_c.pr2$uci<-m3_c.pr$fit+(m3_c.pr2$se*1.96)

# Abund 25 fire + location predictions

# For abundance model predictions (abund_25 and abund_5), predict on link scale and back-transform with exp(); neg. binomial models have log link function. For all other models, use the response scale. 

summary(m3_b)

locationfire.pr<-data.frame(location = rep(factor(levels(sum_dat$location),levels = levels(sum_dat$location)),3),fire_cat = c(rep(levels(sum_dat$fire_cat)[1],2),rep(levels(sum_dat$fire_cat)[2],2),rep(levels(sum_dat$fire_cat)[3],2)))

# fire + location predictions for models to plot graphs, use m3_b

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

m4.set<-list("fire x location"= m4_a, "location"= m4_b, "fire"= m4_c, "null"= m4_d)
m4.tab<-aictab(cand.set = m4.set, second.ord = T, sort = T)

#write.table(m4.tab,file="m4_tab.txt", quote = F, sep = "\t", row.names=F)

# abund_5 predictions for models to plot graphs
fireonly.pr<-data.frame(fire_cat = factor(levels(sum_dat$fire_cat),levels = levels(sum_dat$fire_cat)))
m4_c.pr<-predict(object = m4_c, newdata = fireonly.pr,type = "link", se.fit = T)
m4_c.pr2<-data.frame(fireonly.pr)
m4_c.pr2$fit<-m4_c.pr$fit
m4_c.pr2$se<-m4_c.pr$se
m4_c.pr2$lci<-m4_c.pr$fit-(m4_c.pr2$se*1.96)
m4_c.pr2$uci<-m4_c.pr$fit+(m4_c.pr2$se*1.96)

m4_c.pr2$fit<-exp(m4_c.pr2$fit)
m4_c.pr2$se<-exp(m4_c.pr2$se)
m4_c.pr2$lci<-exp(m4_c.pr2$lci)
m4_c.pr2$uci<-exp(m4_c.pr2$uci)


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

m5.set<-list("fire x location"= m5_a, "location"= m5_b, "fire"= m5_c, "null"= m5_d)
m5.tab<-aictab(cand.set = m5.set, second.ord = T, sort = T)

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

m5_b.pr<-predict(object = m5_b, newdata = locationfire.pr,type = "response", se.fit = T)
m5_b.pr2<-data.frame(locationfire.pr)
m5_b.pr2$fit<-m5_b.pr$fit
m5_b.pr2$se<-m5_b.pr$se
m5_b.pr2$lci<-m5_b.pr$fit-(m5_b.pr2$se*1.96)
m5_b.pr2$uci<-m5_b.pr$fit+(m5_b.pr2$se*1.96)

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

m6.set<-list("fire x location"= m6_a, "location"= m6_b, "fire"= m6_c, "null"= m6_d)
m6.tab<-aictab(cand.set = m6.set, second.ord = T, sort = T)

# species richness of RARE SPECIES, as a function of fire category (generalized linear model with negative binomial linear structure)
head(sum_dat,3); dim(sum_dat)

# Rare species Lowest 25 %
m7_a <- glm.nb(sr_25~fire_cat*location,data = sum_dat)
m7_b <- glm.nb(sr_25~fire_cat+location,data = sum_dat)
m7_c <- glm.nb(sr_25~fire_cat,data = sum_dat)
m7_d <- glm.nb(sr_25~1,data = sum_dat)

#Likelihood ratio test - compares more complicated model (interaction model) with simpler model 
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

m7_c.pr<-predict(object = m7_c, newdata = fireonly.pr,type = "response", se.fit = T)
m7_c.pr2<-data.frame(fireonly.pr)
m7_c.pr2$fit<-m7_c.pr$fit
m7_c.pr2$se<-m7_c.pr$se
m7_c.pr2$lci<-m7_c.pr$fit-(m7_c.pr2$se*1.96)
m7_c.pr2$uci<-m7_c.pr$fit+(m7_c.pr2$se*1.96)

# Rare species 5 %
m8_a <- glm.nb(sr_5~fire_cat*location,data = sum_dat)
m8_b <- glm.nb(sr_5~fire_cat+location,data = sum_dat)
m8_c <- glm.nb(sr_5~fire_cat,data = sum_dat)
m8_d <- glm.nb(sr_5~1,data = sum_dat)

#Likelihood ratio test - compares more complicated model (interaction model) with simpler model 
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
fireonly.pr

m8_c.pr<-predict(object = m8_c, newdata = fireonly.pr,type = "response", se.fit = T)
m8_c.pr2<-data.frame(fireonly.pr)
m8_c.pr2$fit<-m8_c.pr$fit
m8_c.pr2$se<-m8_c.pr$se
m8_c.pr2$lci<-m8_c.pr$fit-(m8_c.pr2$se*1.96)
m8_c.pr2$uci<-m8_c.pr$fit+(m8_c.pr2$se*1.96)

# make AICc table of for all response variables:

m1.tab2<-data.frame(response="Species Richness",m1.tab)
m2.tab2<-data.frame(response="Simpson's Index",m2.tab)
m3.tab2<-data.frame(response="Abund_25",m3.tab)
m4.tab2<-data.frame(response="Abund_5",m4.tab)
m5.tab2<-data.frame(response="Shannon's Index",m5.tab)
m6.tab2<-data.frame(response="Evenness",m6.tab)
m7.tab2<-data.frame(response="Richness Lowest 25%",m7.tab)
m8.tab2<-data.frame(response="Richness Max. 5%",m8.tab)

# Haven't added m7 and m8 to combi.tab yet (richness of rare species):
combi.tab<-rbind(m1.tab2, m2.tab2, m3.tab2, m4.tab2, m5.tab2, m6.tab2, m7.tab2, m8.tab2)

#write.table(combi.tab,file="megatable.txt", quote = F, sep = "\t", row.names=F)

#save.image("04_workspaces/analysed_data.RData")


#Contrasts

#There are 3 contrasts for three categories (u:m, u:b, m:b)

summary(m1_c) #glm.nb model

# fireonly.pr<-data.frame(fire_cat = factor(levels(sum_dat$fire_cat),levels = levels(sum_dat$fire_cat)))

fireonly_contrast<-data.frame(contrast=paste(combn(fireonly.pr$fire_cat,2)[1,],combn(fireonly.pr$fire_cat,2)[2,],sep=':'))

#Create unique model matrix

m1c_conmod <- lm(formula = sp_rich ~ fire_cat, data = sum_dat,x=T)$x

umm_m1c <- unique(m1c_conmod)

rownames(umm_m1c) <- 1:nrow(umm_m1c)

#WARNING NUMERIC SUBSETS - put them in the natural order 2016-2019, c to r

umm_m1c <- umm_m1c[c(1,3,2),]

rownames(umm_m1c) <- 1:nrow(umm_m1c)

#Create a difference matrix

#Each row must be a vector with a length equal to the number of rows in the unique model matrix (umm), e.g. four rows in umm_form matrix will give 6 contrasts. Each row will specify one contrast.

diffm_col <- rbind(
  
  c(-1,1,0,0),
  
  c(-1,0,1,0),
  
  c(-1,0,0,1),
  
  c(0,-1,1,0),
  
  c(0,-1,0,1),
  
  c(0,0,-1,1)
  
)



#Now we have a unique model matrix

umm_col2



#and we have a difference matrix

diffm_col



#and we have the names for the contrast

col_contrast



#calculate the differences and CI's (abun)

col_diff<-data.frame(contrast=col_contrast,diff.est(model = Colinvdiv_mod1,unique.mod.mat = umm_col2,diff.matrix = diffm_col))

col_diff$diff <- ifelse(sign(col_diff$lci)==sign(col_diff$uci),1,0)





#### RARE SPECIES assemblage
#### Species richness of rare species combined with abundance of rare species:

# Lowest 25% (14 species). This metric classified species as rare if their abundances fell into the lowest 25% of the assemblage. For these 14 species, we modelled both the richness of this assemblage and the total abundance of all 14 species. 

# 5 % of max abundance ()
m7_c.pr2
m8_c.pr2

dev.new(width=7, height=6, dpi=80, pointsize=14, noRStudioGD = T)
par(mfrow=c(2,2), mar=c(4.5,4,1,1), mgp=c(2.8,0.8,0), oma=c(0,0,1,6))

# species richness (5% of max.)
plot(c(1:3),m8_c.pr2$fit, xlim=c(0.5,3.5), pch=20, xaxt="n",ylim= c((min(m8_c.pr2$lci)),max(m8_c.pr2$uci)),ylab="Richness (5% of max.)",xlab="", las = 1, cex = 2.5)
arrows(c(1:3),m8_c.pr2$lci,c(1:3),m8_c.pr2$uci,length=0.03,code=3,angle=90)
axis(1,at=c(1:3),labels=F)
axis(1,at=c(0.8,2,3.2), cex.axis=1,labels=m8_c.pr2$fire_cat,tick=F)
title(mgp=c(2.3,0.8,0),xlab="Fire Category")
m8.tab2
text(2.1,max(m8_c.pr2$uci),as.expression(bquote(Delta~"AICc ="~.(paste(round(m8.tab2$Delta_AICc[m8.tab2$Modnames=="fire"]-m1.tab2$Delta_AICc[m8.tab2$Modnames=="null"],2),sep="")))),adj=0,col="red", cex=0.9)
mtext(text="(a)", side = 3, line = 0.5, adj = 0, cex = 1)

# species richness (Lowest 25%)
plot(c(1:3),m7_c.pr2$fit, xlim=c(0.5,3.5), pch=20, xaxt="n",ylim= c((min(m7_c.pr2$lci)),max(m7_c.pr2$uci)),ylab="Richness (lowest 25%)",xlab="", las = 1, cex = 2.5)
arrows(c(1:3),m7_c.pr2$lci,c(1:3),m7_c.pr2$uci,length=0.03,code=3,angle=90)
axis(1,at=c(1:3),labels=F)
axis(1,at=c(0.8,2,3.2),labels=m7_c.pr2$fire_cat,tick=F, cex.axis=1)
title(mgp=c(2.3,0.8,0),xlab="Fire Category")
m7.tab2
text(2.1,max(m7_c.pr2$uci),as.expression(bquote(Delta~"AICc ="~.(paste(round(m7.tab2$Delta_AICc[m7.tab2$Modnames=="fire"]-m7.tab2$Delta_AICc[m7.tab2$Modnames=="null"],2),sep="")))),adj=0,col="red", cex=0.9)
mtext(text="(b)", side = 3, line = 0.5, adj = 0, cex = 1)

par(xpd=NA)
legend(x=4,y=max(m7_c.pr2$uci)+0.1, title = "Sites", legend = c("Fire only", "Hincks","Pinks"), pt.cex = 1.5, pch = c(16, 15, 17), bty = "n", title.adj=0)
par(xpd=F)

# Abundance (5% of max.)
plot(c(1:3),m4_c.pr2$fit, xlim=c(0.5,3.5), pch=20, xaxt="n",ylim= c((min(m4_c.pr2$lci)),max(m4_c.pr2$uci)),ylab="Abundance (5% of max.)",xlab="", las = 1, cex = 2.5)
arrows(c(1:3),m4_c.pr2$lci,c(1:3),m4_c.pr2$uci,length=0.03,code=3,angle=90)
axis(1,at=c(1:3),labels=F)
axis(1,at=c(0.8,2,3.2),cex.axis=1,labels=m4_c.pr2$fire_cat,tick=F)
title(mgp=c(2.3,0.8,0),xlab="Fire Category")
text(2.1,max(m4_c.pr2$uci),as.expression(bquote(Delta~"AICc ="~.(paste(round(m4.tab2$Delta_AICc[m4.tab2$Modnames=="fire"]-m4.tab2$Delta_AICc[m4.tab2$Modnames=="null"],2),sep="")))),adj=0,col="red", cex=0.9)
mtext(text="(c)", side = 3, line = 0.5, adj = 0, cex = 1)

# Abundance (Lowest 25%)
plot(c(1:3)-0.15,m3_b.pr2$fit[m3_b.pr2$location=="Hincks"], xlim=c(0.5,3.5), pch=15, xaxt="n",ylim= c((min(m3_b.pr2$lci)),max(m3_b.pr2$uci)),ylab="Abundance (lowest 25%)",xlab="", las = 1, cex = 1.5)
points(c(1:3)+0.15,m3_b.pr2$fit[m3_b.pr2$location=="Pinks"], xlim=c(0.5,3.5), pch=17, cex = 1.5)
arrows(c(1:3)-0.15,m3_b.pr2$lci[m3_b.pr2$location=="Hincks"],c(1:3)-0.15,m3_b.pr2$uci[m3_b.pr2$location=="Hincks"],length=0.03,code=3,angle=90)
arrows(c(1:3)+0.15,m3_b.pr2$lci[m3_b.pr2$location=="Pinks"],c(1:3)+0.15,m3_b.pr2$uci[m3_b.pr2$location=="Pinks"],length=0.03,code=3,angle=90)
axis(1,at=c(1:3),labels=F)
axis(1,at=c(0.8,2,3.2), cex.axis=1,labels=m3_b.pr2$fire_cat[m3_b.pr2$location=="Hincks"],tick=F)
title(mgp=c(2.3,0.8,0),xlab="Fire Category")
text(2.1,max(m3_b.pr2$uci),as.expression(bquote(Delta~"AICc ="~.(paste(round(m3.tab2$Delta_AICc[m3.tab2$Modnames=="fire"]-m3.tab2$Delta_AICc[m3.tab2$Modnames=="null"],2),sep="")))),adj=0,col="red", cex=0.9)
mtext(text="(d)", side = 3, line = 0.5, adj = 0, cex = 1)


#### ALL SPECIES assemblage
#### Species richness all species:

dev.new(width=7, height=6, dpi=80, pointsize=14, noRStudioGD = T)
par(mfrow=c(2,2), mar=c(4.5,4,1,1), mgp=c(2.8,0.8,0), oma=c(0,0,1,6))

# species richness
plot(c(1:3),m1_c.pr2$fit, xlim=c(0.5,3.5), pch=20, xaxt="n",ylim= c((min(m1_c.pr2$lci)),max(m1_c.pr2$uci)),ylab="Species Richness",xlab="", las = 1, cex = 2.5)
arrows(c(1:3),m1_c.pr2$lci,c(1:3),m1_c.pr2$uci,length=0.03,code=3,angle=90)
axis(1,at=c(1:3),labels=F)
axis(1,at=c(0.8,2,3.2), cex.axis=1,labels=m1_c.pr2$fire_cat,tick=F)
title(mgp=c(2.3,0.8,0),xlab="Fire Category")
m1.tab2
text(2.1,max(m1_c.pr2$uci),as.expression(bquote(Delta~"AICc ="~.(paste(round(m1.tab2$Delta_AICc[m1.tab2$Modnames=="fire"]-m1.tab2$Delta_AICc[m1.tab2$Modnames=="null"],2),sep="")))),adj=0,col="red",cex=0.9)
mtext(text="(a)", side = 3, line = 0.5, adj = 0, cex = 1)

# simps diversity index, no changes required
plot(c(1:3),m2_c.pr2$fit, xlim=c(0.5,3.5), pch=20, xaxt="n",ylim= c((min(m2_c.pr2$lci)),max(m2_c.pr2$uci)),ylab="Simpson's Diversity Index",xlab="", las = 1, cex = 2.5)
arrows(c(1:3),m2_c.pr2$lci,c(1:3),m2_c.pr2$uci,length=0.03,code=3,angle=90)
axis(1,at=c(1:3),labels=F)
axis(1,at=c(0.8,2,3.2), cex.axis=1,labels=m2_c.pr2$fire_cat,tick=F)
title(mgp=c(2.3,0.8,0),xlab="Fire Category")
text(2.1,max(m2_c.pr2$uci),as.expression(bquote(Delta~"AICc ="~.(paste(round(m2.tab2$Delta_AICc[m2.tab2$Modnames=="fire"]-m2.tab2$Delta_AICc[m2.tab2$Modnames=="null"],2),sep="")))),adj=0,col="red", cex=0.9)
mtext(text="(b)", side = 3, line = 0.5, adj = 0, cex = 1)

par(xpd=NA)
legend(x=4,y=max(m2_c.pr2$uci)+0.1, title = "Sites", legend = c("Fire only", "Hincks","Pinks"), pt.cex = 1.5, pch = c(16, 15, 17), bty = "n", title.adj=0)
par(xpd=F)

# shann_ind plots for fire+location

plot(c(1:3)-0.15,m5_b.pr2$fit[m5_b.pr2$location=="Hincks"], xlim=c(0.5,3.5), pch=15, xaxt="n",ylim= c((min(m5_b.pr2$lci)),max(m5_b.pr2$uci)+0.03),ylab="Shannon's Index",xlab="", las = 1, cex = 1.5)
points(c(1:3)+0.15,m5_b.pr2$fit[m5_b.pr2$location=="Pinks"], xlim=c(0.5,3.5), pch=17, cex = 1.5)
arrows(c(1:3)-0.15,m5_b.pr2$lci[m5_b.pr2$location=="Hincks"],c(1:3)-0.15,m5_b.pr2$uci[m5_b.pr2$location=="Hincks"],length=0.03,code=3,angle=90)
arrows(c(1:3)+0.15,m5_b.pr2$lci[m5_b.pr2$location=="Pinks"],c(1:3)+0.15,m5_b.pr2$uci[m5_b.pr2$location=="Pinks"],length=0.03,code=3,angle=90)
axis(1,at=c(1:3),labels=F)
axis(1,at=c(0.8,2,3.2),labels=m5_b.pr2$fire_cat[m5_b.pr2$location=="Hincks"],tick=F)
title(mgp=c(2.3,0.8,0),xlab="Fire Category")
text(2,max(m5_b.pr2$uci)+0.03,as.expression(bquote(Delta~"AICc ="~.(paste(round(m5.tab2$Delta_AICc[m5.tab2$Modnames=="fire"]-m5.tab2$Delta_AICc[m5.tab2$Modnames=="null"],2),sep="")))),adj=0,col="dark green", cex=0.9)
mtext(text="(c)", side = 3, line = 0.5, adj = 0, cex = 1)

# evenness plots for fire
plot(c(1:3),m6_c.pr2$fit, xlim=c(0.5,3.5), pch=20, xaxt="n",ylim= c((min(m6_c.pr2$lci)),max(m6_c.pr2$uci)),ylab="Evenness",xlab="", las = 1, cex = 2.5)
arrows(c(1:3),m6_c.pr2$lci,c(1:3),m6_c.pr2$uci,length=0.03,code=3,angle=90)
axis(1,at=c(1:3),labels=F)
axis(1,at=c(0.8,2,3.2),labels=m6_c.pr2$fire_cat,tick=F)
title(mgp=c(2.3,0.8,0),xlab="Fire category")
text(2,max(m6_c.pr2$uci),as.expression(bquote(Delta~"AICc ="~.(paste(round(m6.tab2$Delta_AICc[m6.tab2$Modnames=="fire"]-m6.tab2$Delta_AICc[m6.tab2$Modnames=="null"],2),sep="")))),adj=0,col="dark green", cex=0.9)
mtext(text="(d)", side = 3, line = 0.5, adj = 0, cex = 1)







