#### LOCATION ONLY plot code

#data frame for location only predictions

locationonly.pr<-data.frame(location = factor(levels(sum_dat$location),levels = levels(sum_dat$location)))

# location predictions for models to plot graphs, use m1_b

m1_b.pr<-predict(object = m1_b, newdata = locationonly.pr,type = "response", se.fit = T)
m1_b.pr2<-data.frame(locationonly.pr)
m1_b.pr2$fit<-m1_b.pr$fit
m1_b.pr2$se<-m1_b.pr$se
m1_b.pr2$lci<-m1_b.pr$fit-(m1_b.pr2$se*1.96)
m1_b.pr2$uci<-m1_b.pr$fit+(m1_b.pr2$se*1.96)

#data frame for location only predictions

locationonly.pr<-data.frame(location = factor(levels(sum_dat$location),levels = levels(sum_dat$location)))

# location predictions for models to plot graphs, use m2_b

m2_b.pr<-predict(object = m2_b, newdata = locationonly.pr,type = "response", se.fit = T)
m2_b.pr2<-data.frame(locationonly.pr)
m2_b.pr2$fit<-m2_b.pr$fit
m2_b.pr2$se<-m2_b.pr$se
m2_b.pr2$lci<-m2_b.pr$fit-(m2_b.pr2$se*1.96)
m2_b.pr2$uci<-m2_b.pr$fit+(m2_b.pr2$se*1.96)

#data frame for location only predictions

locationonly.pr<-data.frame(location = factor(levels(sum_dat$location),levels = levels(sum_dat$location)))

# location predictions for models to plot graphs, use m2_b

m4_b.pr<-predict(object = m4_b, newdata = locationonly.pr,type = "response", se.fit = T)
m4_b.pr2<-data.frame(locationonly.pr)
m4_b.pr2$fit<-m4_b.pr$fit
m4_b.pr2$se<-m4_b.pr$se
m4_b.pr2$lci<-m4_b.pr$fit-(m4_b.pr2$se*1.96)
m4_b.pr2$uci<-m4_b.pr$fit+(m4_b.pr2$se*1.96)




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

# plots for strongest fit

dev.new(width=7, height=9, dpi=75, pointsize=16, noRStudioGD = T)
par(mfrow=c(3,2), mar=c(4.5,4,1,1), mgp=c(2.8,0.8,0), oma=c(0,0,1,6))

# species richness
plot(c(1:3),m1_c.pr2$fit, xlim=c(0.5,3.5), pch=20, xaxt="n",ylim= c((min(m1_c.pr2$lci)),max(m1_c.pr2$uci)),ylab="Species Richness",xlab="", las = 1, cex = 2.5)
arrows(c(1:3),m1_c.pr2$lci,c(1:3),m1_c.pr2$uci,length=0.03,code=3,angle=90)
axis(1,at=c(1:3),labels=F)
axis(1,at=c(1,2,3),labels=m1_c.pr2$fire_cat,tick=F)
title(mgp=c(2.3,0.8,0),xlab="Fire Category")
m1.tab2
text(2.1,max(m1_c.pr2$uci),as.expression(bquote(Delta~"AICc ="~.(paste(round(m1.tab2$Delta_AICc[m1.tab2$Modnames=="fire"]-m1.tab2$Delta_AICc[m1.tab2$Modnames=="null"],2),sep="")))),adj=0,col="red")
mtext(text="(a)", side = 3, line = 0.5, adj = 0, cex = 1)

# simps diversity index, no changes required
plot(c(1:3),m2_c.pr2$fit, xlim=c(0.5,3.5), pch=20, xaxt="n",ylim= c((min(m2_c.pr2$lci)),max(m2_c.pr2$uci)),ylab="Simpson's Diversity Index",xlab="", las = 1, cex = 2.5)
arrows(c(1:3),m2_c.pr2$lci,c(1:3),m2_c.pr2$uci,length=0.03,code=3,angle=90)
axis(1,at=c(1:3),labels=F)
axis(1,at=c(1,2,3),labels=m2_c.pr2$fire_cat,tick=F)
title(mgp=c(2.3,0.8,0),xlab="Fire Category")
text(2.1,max(m2_c.pr2$uci),as.expression(bquote(Delta~"AICc ="~.(paste(round(m2.tab2$Delta_AICc[m2.tab2$Modnames=="fire"]-m2.tab2$Delta_AICc[m2.tab2$Modnames=="null"],2),sep="")))),adj=0,col="red")
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
axis(1,at=c(1,2,3),labels=m5_b.pr2$fire_cat[m5_b.pr2$location=="Hincks"],tick=F)
title(mgp=c(2.3,0.8,0),xlab="Fire Category")
text(2,max(m5_b.pr2$uci)+0.03,as.expression(bquote(Delta~"AICc ="~.(paste(round(m5.tab2$Delta_AICc[m5.tab2$Modnames=="fire"]-m5.tab2$Delta_AICc[m5.tab2$Modnames=="null"],2),sep="")))),adj=0,col="dark green")
mtext(text="(c)", side = 3, line = 0.5, adj = 0, cex = 1)

# evenness plots for fire
plot(c(1:3),m6_c.pr2$fit, xlim=c(0.5,3.5), pch=20, xaxt="n",ylim= c((min(m6_c.pr2$lci)),max(m6_c.pr2$uci)),ylab="Evenness",xlab="", las = 1, cex = 2.5)
arrows(c(1:3),m6_c.pr2$lci,c(1:3),m6_c.pr2$uci,length=0.03,code=3,angle=90)
axis(1,at=c(1:3),labels=F)
axis(1,at=c(1,2,3),labels=m6_c.pr2$fire_cat,tick=F)
title(mgp=c(2.3,0.8,0),xlab="Fire category")
text(2,max(m6_c.pr2$uci),as.expression(bquote(Delta~"AICc ="~.(paste(round(m6.tab2$Delta_AICc[m6.tab2$Modnames=="fire"]-m6.tab2$Delta_AICc[m6.tab2$Modnames=="null"],2),sep="")))),adj=0,col="dark green")
mtext(text="(d)", side = 3, line = 0.5, adj = 0, cex = 1)

# abund_5 plot
plot(c(1:3),m4_c.pr2$fit, xlim=c(0.5,3.5), pch=20, xaxt="n",ylim= c((min(m4_c.pr2$lci)),max(m4_c.pr2$uci)),ylab="5% of max abundance",xlab="", las = 1, cex = 2.5)
arrows(c(1:3),m4_c.pr2$lci,c(1:3),m4_c.pr2$uci,length=0.03,code=3,angle=90)
axis(1,at=c(1:3),labels=F)
axis(1,at=c(1,2,3),labels=m4_c.pr2$fire_cat,tick=F)
title(mgp=c(2.3,0.8,0),xlab="Fire Category")
text(2.1,max(m4_c.pr2$uci),as.expression(bquote(Delta~"AICc ="~.(paste(round(m4.tab2$Delta_AICc[m4.tab2$Modnames=="fire"]-m4.tab2$Delta_AICc[m4.tab2$Modnames=="null"],2),sep="")))),adj=0,col="red")
mtext(text="(e)", side = 3, line = 0.5, adj = 0, cex = 1)

# abund_25 plot fire*location

plot(c(1:3)-0.15,m3_b.pr2$fit[m3_b.pr2$location=="Hincks"], xlim=c(0.5,3.5), pch=15, xaxt="n",ylim= c((min(m3_b.pr2$lci)),max(m3_b.pr2$uci)),ylab="Lowest 25%",xlab="", las = 1, cex = 1.5)
points(c(1:3)+0.15,m3_b.pr2$fit[m3_b.pr2$location=="Pinks"], xlim=c(0.5,3.5), pch=17, cex = 1.5)
arrows(c(1:3)-0.15,m3_b.pr2$lci[m3_b.pr2$location=="Hincks"],c(1:3)-0.15,m3_b.pr2$uci[m3_b.pr2$location=="Hincks"],length=0.03,code=3,angle=90)
arrows(c(1:3)+0.15,m3_b.pr2$lci[m3_b.pr2$location=="Pinks"],c(1:3)+0.15,m3_b.pr2$uci[m3_b.pr2$location=="Pinks"],length=0.03,code=3,angle=90)
axis(1,at=c(1:3),labels=F)
axis(1,at=c(1,2,3),labels=m3_b.pr2$fire_cat[m3_b.pr2$location=="Hincks"],tick=F)
title(mgp=c(2.3,0.8,0),xlab="Fire Category")
text(2.1,max(m3_b.pr2$uci),as.expression(bquote(Delta~"AICc ="~.(paste(round(m3.tab2$Delta_AICc[m3.tab2$Modnames=="fire"]-m3.tab2$Delta_AICc[m3.tab2$Modnames=="null"],2),sep="")))),adj=0,col="red")
mtext(text="(f)", side = 3, line = 0.5, adj = 0, cex = 1)