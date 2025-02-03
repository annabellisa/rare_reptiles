# ------------------------------------ #
# ---------- RARE REPTILES  ---------- #
# ------------------------------------ #

### Analysis fire effects on rarity and dominance in a reptile community 

# ---------- 05_range_size  ---------- #

### Script authors: Annabel Smith

# If starting fresh:
# load("04_workspaces/analysed_data.RData"); rm(list = setdiff(ls(), c("sp_sum2")))

# Load workspace: 
load("04_workspaces/range_size.RData")

# save.image("04_workspaces/range_size.RData")

# Load functions:
invisible(lapply(paste("02_functions/",dir("02_functions"),sep=""), function(x) source(x)))

# Extract geographical range size data:

# Meiri S (2024). SquamBase—A database of squamate (Reptilia: Squamata) traits. Global Ecology and Biogeography 33, e13812.

# SquamBase:

# Should be 11744 rows (excl. header)
# 84 columns

dir("01_data")

sbase<-read.table("01_data/Meiri_2024_GEB_squamBase1.txt",header=T, strip.white = T,quote ="",sep="\t",check.names = F, fill=T)

# species data 
head(sp_sum2);dim(sp_sum2)

sp.all<-sp_sum2$Species_Name
sp.geo<-data.frame(Species=sp_sum2$Species_Name)

table(sp.all %in% sbase$"Species name (Binomial)")

# Extract data for species with matching names:

# Both measures are in km^2
geo<-data.frame(species=sp.all[sp.all %in% sbase$"Species name (Binomial)"],sbase[which(sbase$"Species name (Binomial)" %in% sp.all),suppressWarnings(grep("Geographic",colnames(sbase)))])
colnames(geo)<-c("Species","geoCaetano","geoIUCN")
head(geo)

# Merge with full dataset:

sp.geo<-merge(sp.geo,geo,by="Species",all.x=T, all.y=F)
head(sp.geo); dim(sp.geo)

# Find species without matching names:
sp.all[!sp.all %in% sbase$"Species name (Binomial)"]
no.match<-sp.all[!sp.all %in% sbase$"Species name (Binomial)"]
no.match<-no.match[!no.match=="Pogona spp"]

# Parasuta spectabilis = Suta spectabilis

"Suta spectabilis" %in% sbase$"Species name (Binomial)"

sp.geo[which(sp.geo$Species=="Parasuta spectabilis"),c("geoCaetano","geoIUCN")]<-sbase[which(sbase$"Species name (Binomial)" %in% "Suta spectabilis"),suppressWarnings(grep("Geographic",colnames(sbase)))]

# Tympanocryptis lineata ssp. lineata = Tympanocryptis lineata

"Tympanocryptis lineata" %in% sbase$"Species name (Binomial)"

sp.geo[which(sp.geo$Species=="Tympanocryptis lineata ssp. lineata"),c("geoCaetano","geoIUCN")]<-sbase[which(sbase$"Species name (Binomial)" %in% "Tympanocryptis lineata"),suppressWarnings(grep("Geographic",colnames(sbase)))]

# Ramphotyphlops bicolor = Anilios bicolor 

"Anilios bicolor" %in% sbase$"Species name (Binomial)"

sp.geo[which(sp.geo$Species=="Ramphotyphlops bicolor"),c("geoCaetano","geoIUCN")]<-sbase[which(sbase$"Species name (Binomial)" %in% "Anilios bicolor"),suppressWarnings(grep("Geographic",colnames(sbase)))]

# Ramphotyphlops bituberculatus = Anilios bituberculatus 

"Anilios bituberculatus" %in% sbase$"Species name (Binomial)"

sp.geo[which(sp.geo$Species=="Ramphotyphlops bituberculatus"),c("geoCaetano","geoIUCN")]<-sbase[which(sbase$"Species name (Binomial)" %in% "Anilios bituberculatus"),suppressWarnings(grep("Geographic",colnames(sbase)))]

head(sp.geo)

sp.geo$newname<-sp.geo$Species
sp.geo$newname[which(sp.geo$Species %in% no.match)]<-c("Suta spectabilis","Anilios bicolor","Anilios bituberculatus","Tympanocryptis lineata")

# save.image("04_workspaces/range_size.RData")

# Get Pogona range sizes
# P. barbata, P. vittatus, P. minor

sbase[which(sbase$"Species name (Binomial)" %in% c("Pogona barbata", "Pogona vitticeps", "Pogona minor")),c(1,suppressWarnings(grep("Geographic",colnames(sbase))))]

sp.geo[which(sp.geo$Species=="Pogona spp"),c(2,3)]<-">1600000"

# Get other potentially useful info:

trait<-data.frame(Species=sp.geo$newname)
trait<-trait[-which(trait$Species=="Pogona spp"),]
trait<-data.frame(Species=c(trait,"Pogona barbata", "Pogona vitticeps", "Pogona minor"))
trait<-data.frame(Species=trait[order(trait$Species),])

table(trait$Species %in% sbase$"Species name (Binomial)")

head(sbase[,1:5],2); dim(sbase)

trait.all<-sbase[which(sbase$"Species name (Binomial)" %in% trait$Species),]
head(trait.all[,1:5],2); dim(trait.all)

# write.table(trait.all,file="trait.all.txt",quote=F,sep="\t",row.names=F)
# write.table(sp.geo,file="sp.geo.txt",quote=F,sep="\t",row.names=F)

# save.image("04_workspaces/range_size.RData")




