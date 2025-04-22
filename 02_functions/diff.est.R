# This function estimates differences, SEs and CIs among all pairs of levels within factorial variables (least-squares mean contrasts) on the link scale. 

# I wrote it for logistic and negbin models and later adapted it for lmer. In Nov 2021, I adapted it for gam, glmer and glmmadmb. In Sept 2023, I adapted it for negbin models (glm.nb), for Amber Lim's rare_reptiles. 

# Inputs are the model to estimate from, a unique model matrix for the parameters of interest (with other params set to zero, and the differences matrix which specifies how to calculate the differences):

# Author: Annabel Smith & Wade Blanchard

diff.est<-function(model,unique.mod.mat,diff.matrix){
  
  # ESTIMATES:
  
  contrasts<-diff.matrix%*%unique.mod.mat
  
  if(class(model)[1]=="glmerMod") mod.coef<-summary(model)$coefficients[,1]
  
  if(class(model)[1]=="glmmadmb") mod.coef<-summary(model)$coefficients[,1]
  
  if(class(model)[1]=="gam") mod.coef<-summary(model)$p.table[,1]
  
  if(class(model)[1]=="negbin") mod.coef<-summary(model)$coefficients[,1]
  
  if(class(model)[1]=="glm") mod.coef<-summary(model)$coefficients[,1]
  
  diffs<-contrasts%*%mod.coef
  
  # SE:
  
  mod.vcov<-vcov(model)
  
  diff.se<-sqrt(diag(contrasts%*%mod.vcov%*%t(contrasts)))
  
  # CI
  
  ndf<-data.frame(est=diffs, se=diff.se)
  
  ndf$lci<-ndf$est-(ndf$se*1.96)
  
  ndf$uci<-ndf$est+(ndf$se*1.96)
  
  return(ndf)
  
}