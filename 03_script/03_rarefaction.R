
# 03_rarefaction

# Load and tidy workspace and remove everything except necessary objects:
load("04_workspaces/analysed_data.RData"); rm(list=setdiff(ls(), c("sum_dat","rare.data","sp_div2","sp_div_25","sp_div_5")))

# Load functions:
invisible(lapply(paste("02_functions/",dir("02_functions"),sep=""), function(x) source(x)))

# save.image("04_workspaces/rarefaction.RData")

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






