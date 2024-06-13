# ################################# ####
# The legacy of Novel Communities   ####
# Author:    Timothy L Staples      ####
# Collaborators: John Pandolfi      #### 
#                Wolfgang Kiessling ####
# ################################# ####
# SET-UP ####
# Global attributes & working directories ####

rm(list=ls())

setwd("/Users/uqtstapl/Library/CloudStorage/Dropbox/Tim/Post-doc/Research projects/novel_comms_legacy/prodCode")

# Packages & functions ####

# source functions from 'functions' sub-folder
sapply(paste0("./functions/", list.files("./functions", pattern =".R")), source)

package.loader(c("vegan", "DHARMa", "lme4", "multcomp", "plotrix", "shape", "hilldiv",
                 "merTools", "performance", "abind", "sf", "rworldmap", "readxl",
                 "terra", "ncdf4", "tidyr", "raster", "performance"))

# A little custom function to add dates to output files
date.wrap <- function(string, ext){
  paste0(string, " ", Sys.Date(), ext)
}

# converts x and y values into angles on a continuous 360 degree scale
cont.angles <- function(x,y, offset = 0){
  angles <- -atan2(x, y) * (180/pi)
  angles <- ifelse(angles > 0, angles, 360 - abs(angles))
  rot.angles <- angles + 90 + offset
  rot.angles <- ifelse(rot.angles >= 360, rot.angles - 360, rot.angles)
  return(rot.angles)
}

# function to lighten or darken colours using ramps
colorShader <- function(col, scale, direction = "lighten"){
  rgb(colorRamp(c(col, ifelse(direction=="lighten", "white", "black")))(scale)/255)
}

# convert data to scale (default 0-1). For colour ramps for plotting.
unitScale <- function(x, custMin = NULL, custMax=NULL){
  if(is.null(custMax)){custMax = max(x, na.rm=TRUE)}
  if(is.null(custMin)){custMin = min(x, na.rm=TRUE)}
  (x - custMin) / (custMax - custMin)
}
# ------------------ ####
# INITIAL PROCESSING ####

plant.record.df <- readRDS("./rawdata/processedRecords.rds")
plant.record.df <- droplevels(plant.record.df[plant.record.df$elementtype == "pollen",])

# ------------------ ####
#           Get time-series location ####

world <- getMap("high")

site.df <- plant.record.df[!duplicated(plant.record.df$siteid) & complete.cases(plant.record.df[,c("long","lat")]),]
#coordinates(site.coords) <- c("long", "lat")
site.coords.raw <- site.df[,c("long", "lat")]
coordinates(site.coords.raw) <- c("long", "lat")
proj4string(site.coords.raw) <- proj4string(world)

site.df.save <- site.df
site.df <- cbind(site.df, site.coords.raw %over% world)

# some of these are missing info, likely because their long/lat are
# slightly off the coast
# we can pick up the nearest polygon and give them those attributes,
# assuming the distance is smaller than some error margin

missing.sites <- site.df$siteid[is.na(site.df$Stern)]
missing.site.df <- site.df.save[site.df.save$siteid %in% missing.sites,]
missing.site.df <- missing.site.df[match(missing.site.df$siteid, missing.sites),]
world = st_as_sf(world)

missing.site.geo <- do.call("rbind", lapply(missing.sites, function(x){
  
  print(x) 
  temp.coords <- site.df[site.df$siteid == x, c("long","lat")] 
  coordinates(temp.coords) = c("long", "lat")
  temp.coords = st_as_sf(temp.coords)
  st_crs(temp.coords) <- st_crs(world)
  
  temp.dist <- st_distance(temp.coords, world)
  min.dist <- which.min(temp.dist)

  good.data <- world[min.dist,]
  return(good.data)
  
}))

missing.site.df <- cbind(missing.site.df, missing.site.geo)
missing.site.df <- missing.site.df[,colnames(missing.site.df) != "geometry"]

site.df <- rbind(site.df[!site.df$siteid %in% missing.site.df$siteid,],
                         missing.site.df)

# remove sites not in Europe or Nth America (poor sampling)
site.df <- droplevels(site.df[site.df$REGION %in% c("Europe", "North America"),])

# remove records not in our continental regions, and samples outside of 25000 ybp
plant.record.df <- droplevels(plant.record.df[plant.record.df$siteid %in% site.df$siteid,])
plant.record.df <- droplevels(plant.record.df[plant.record.df$age <= 25100, ])

# Look for duplicate samples
dupeSamps <- paste(plant.record.df$siteid,
                   plant.record.df$sampleid, 
                   plant.record.df$variablename, 
                   plant.record.df$value, sep=".")

plant.record.df <- droplevels(plant.record.df[!duplicated(dupeSamps),])

plant.record.df = merge(plant.record.df, site.df[,c("siteid", "REGION")],
                         by.x="siteid", by.y="siteid", all.x=TRUE, all.y=FALSE, sort=FALSE)

#           Pollen production adjustments ####

# read in PPE
PPE = as.data.frame(read_excel("./rawdata/RPP_Dataset_v2_Table_5_6.xlsx"))
PPE = PPE[,c(1,2,8,12,16)]
colnames(PPE) = c("type","taxon","ppeAm", "ppeEur", "ppeNthHem")

PPE$taxon[PPE$taxon == "Sambucus nigra-type"] = "Sambucus nigra"
PPE$taxon[PPE$taxon == "Asteraceae"] = "Compositae"
PPE$taxon[PPE$taxon == "Fabaceae"] = "Leguminosae"

# pair as best as possible in the following order - Continent + genus, genus, continent + family, family
taxTable = plant.record.df[,c("genus", "family", "REGION")]
taxTable = taxTable[!duplicated(taxTable[,c("genus", "REGION")]),]

eurTable = taxTable[taxTable$REGION == "Europe",]
eurTable = merge(eurTable, PPE[,c("taxon", "ppeEur", "ppeNthHem")], by.x="genus", by.y = "taxon",
                 all.x=TRUE, all.y=FALSE, sort=FALSE)
colnames(eurTable)[-(1:3)] = c("ppeEurGen", "ppeHemGen")

eurTable = merge(eurTable, PPE[,c("taxon", "ppeEur", "ppeNthHem")], by.x="family", by.y = "taxon",
                 all.x=TRUE, all.y=FALSE, sort=FALSE)
eurPPEraw = eurTable

# sample sizes
firstNona = apply(eurTable[,-(1:3)], 1, function(x){which(!is.na(x))[1]})
table(firstNona)

eurPPE = eurTable[,-(1:3)][cbind(1:nrow(eurTable), apply(eurTable[,-(1:3)], 1, function(x){which(!is.na(x))[1]}))]

eurTable = cbind(eurTable[,1:3],
                 PPE = eurPPE)
eurTable[is.na(eurTable$PPE),]

# Nth Am
amTable = taxTable[taxTable$REGION == "North America",]
amTable = merge(amTable, PPE[,c("taxon", "ppeAm", "ppeNthHem")], by.x="genus", by.y = "taxon",
                 all.x=TRUE, all.y=FALSE, sort=FALSE)
colnames(amTable)[-(1:3)] = c("ppeamGen", "ppeHemGen")

amTable = merge(amTable, PPE[,c("taxon", "ppeAm", "ppeNthHem")], by.x="family", by.y = "taxon",
                 all.x=TRUE, all.y=FALSE, sort=FALSE)
amPPEraw = amTable

# sample sizes

# taxa
firstNona = apply(amTable[,-(1:3)], 1, function(x){which(!is.na(x))[1]})

amPPE = amTable[,-(1:3)][cbind(1:nrow(amTable), apply(amTable[,-(1:3)], 1, function(x){which(!is.na(x))[1]}))]

amTable = cbind(amTable[,1:3],
                 PPE = amPPE)
amTable[is.na(amTable$PPE),]

neoPPE = rbind(eurTable, amTable)
neoPPE = neoPPE[!is.na(neoPPE$family),]

# genSub
genPPE = merge(plant.record.df[!is.na(plant.record.df$genus),], 
               neoPPE[,c("genus","REGION","PPE")], by.x=c("genus", "REGION"), by.y=c("genus", "REGION"),
                all.x=TRUE, all.y=FALSE, sort=FALSE)

neoPPEFam = array2DF(tapply(neoPPE$PPE, list(neoPPE$family, neoPPE$REGION), mean, na.rm=TRUE))
colnames(neoPPEFam) = c("family", "REGION", "PPE")

famPPE = merge(plant.record.df[is.na(plant.record.df$genus),], neoPPEFam[,c("family","REGION","PPE")], by.x=c("family", "REGION"), by.y=c("family", "REGION"),
                        all.x=TRUE, all.y=FALSE, sort=FALSE)
genPPE = genPPE[,match(colnames(famPPE), colnames(genPPE))]

plantPPE = rbind(famPPE, genPPE)

# how many samples do we have PPE for?
ppeStats = tapply(plantPPE$value, is.na(plantPPE$PPE), sum, na.rm=TRUE)
ppeStats / sum(ppeStats)

plantPPE$countPPE = plantPPE$value / plantPPE$PPE

plantPPE = droplevels(plantPPE[!is.na(plantPPE$countPPE),])

# sample sizes by count

continentGenus = as.data.frame(tapply(plant.record.df$value, 
                                      list(plant.record.df$genus,
                                           plant.record.df$REGION),
                        sum))

continentFam = with(plant.record.df[is.na(plant.record.df$genus),],
                      as.data.frame(tapply(value, 
                                      list(family, REGION),
                                      sum)))

# how many have continent specific sums?
contGenCount = sum(continentGenus$Europe[rownames(continentGenus) %in% eurPPEraw$genus[!is.na(eurPPEraw$ppeEurGen)]],
                     continentGenus$`North America`[rownames(continentGenus) %in% amPPEraw$genus[!is.na(amPPEraw$ppeamGen)]])
worldGenCount = sum(continentGenus$Europe[rownames(continentGenus) %in% eurPPEraw$genus[!is.na(eurPPEraw$ppeHemGen) & is.na(eurPPEraw$ppeEurGen)]],
                     continentGenus$`North America`[rownames(continentGenus) %in% amPPEraw$genus[!is.na(amPPEraw$ppeHemGen) & is.na(amPPEraw$ppeamGen)]])

contFamCount = sum(continentGenus$Europe[rownames(continentGenus) %in% eurPPEraw$genus[!is.na(eurPPEraw$ppeEur) & is.na(eurPPEraw$ppeEurGen) & is.na(eurPPEraw$ppeHemGen)]],
                   continentGenus$`North America`[rownames(continentGenus) %in% amPPEraw$genus[!is.na(amPPEraw$ppeAm) & is.na(amPPEraw$ppeHemGen) & is.na(amPPEraw$ppeamGen)]])

contFamCount = sum(continentGenus$Europe[rownames(continentGenus) %in% eurPPEraw$genus[!is.na(eurPPEraw$ppeEur) & is.na(eurPPEraw$ppeEurGen) & is.na(eurPPEraw$ppeHemGen)]],
                   continentGenus$`North America`[rownames(continentGenus) %in% amPPEraw$genus[!is.na(amPPEraw$ppeAm) & is.na(amPPEraw$ppeHemGen) & is.na(amPPEraw$ppeamGen)]])

cumsum(c(contGenCount, worldGenCount, contFamCount)) / sum(plant.record.df$value, na.rm=TRUE)

#           Sample binning ####

plantPPE$sampleBin = cut(plantPPE$age, breaks=seq(-100,25100, 200))

plantPPE$sampleBin = rowMeans(cbind(as.numeric( sub("\\((.+),.*", "\\1", plantPPE$sampleBin)),
                                            as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", plantPPE$sampleBin))))

# remove aquatic families
plantPPE <- droplevels(plantPPE[!plantPPE$family %in% 
                                                  c("Potamogetonaceae", "Nymphaeaceae", "Typhaceae", "Cabombaceae", "Alismataceae",
                                                    "Haloragaceae"),])

plantPPE = droplevels(plantPPE[!is.na(plantPPE$value),])

# sample size calculation ####

# how many have family/genus records
sum(!is.na(plantPPE$family))
tapply(plantPPE$value, !is.na(plantPPE$family), sum, na.rm=TRUE) / sum(plantPPE$value, na.rm=TRUE)

sum(!is.na(plantPPE$genus))
tapply(plantPPE$value, !is.na(plantPPE$genus), sum, na.rm=TRUE) / sum(plantPPE$value, na.rm=TRUE)

# how many have PPE ####
sum(!is.na(plantPPE$PPE))
tapply(plantPPE$value, !is.na(plantPPE$PPE), sum, na.rm=TRUE) / sum(plantPPE$value, na.rm=TRUE)

#           final removing ####

# next we need to examine chronology
summary(!is.na(plantPPE$age))
#remove records with no chronology
plantPPE <- droplevels(plantPPE[!is.na(plantPPE$age), ])

plantPPE$site = plantPPE$siteid
plantPPE <- plantPPE[order(plantPPE$site),]

# remove estimates from the future
summary(plantPPE$age > -70)
plantPPE <- plantPPE[plantPPE$age > -70,]

# remove estimates outside our sample period
plantPPE <- plantPPE[plantPPE$age <= 25100,]

# remove sea cores
plantPPE <- plantPPE[plantPPE$elev >= -10,]

# remove any NA rows that have creeped in
plantPPE <- plantPPE[complete.cases(plantPPE[,c("siteid", "sampleid", "variablename", "value")]),]

saveRDS(plantPPE, paste0("./rawdata/processed_family_records.rds"))
write.csv(site.df, "./rawdata/siteDf.csv")

#           Novelty analysis ####

plantPPE <- droplevels(plantPPE)
plantPPE$rawCount = plantPPE$value
plantSlimmed <- plantPPE[complete.cases(plantPPE[,c("family", "site", "age", "countPPE", "rawCount")]),
                                                          c("family", "site", "age", "countPPE", "rawCount")]
colnames(plantSlimmed)[colnames(plantSlimmed) == "countPPE"] = "count"

all.novel <- neotoma.novelty(dataset = plantSlimmed,
                             ssmat.type = "abund",
                             bins = seq(-100,max(plantSlimmed$age), 200),
                             rich.cutoff = c(100, 5000),
                             age.limits = c(-150, 25100),
                             taxon.res = "family",
                             bin.cutoff = 10,
                             taxa.cutoff = 5,
                             novel.alpha = 0.05,
                             novel.metric = "bray",
                             sqrt.mat=TRUE)

# raw counts novelty
plantSlimmedRaw = plantSlimmed
plantSlimmedRaw$count = plantSlimmedRaw$rawCount
all.novel.raw <- neotoma.novelty(dataset = plantSlimmedRaw,
                             ssmat.type = "abund",
                             bins = seq(-100,max(plantSlimmedRaw$age), 200),
                             rich.cutoff = c(100, 5000),
                             age.limits = c(-150, 25100),
                             taxon.res = "family",
                             bin.cutoff = 10,
                             taxa.cutoff = 5,
                             novel.alpha = 0.05,
                             novel.metric = "bray",
                             sqrt.mat=TRUE)

saveRDS(all.novel, "./outputs/all neotoma novelty.rds")
saveRDS(all.novel.raw, "./outputs/all neotoma novelty.raw.rds")

# now do the same thing but with genus-level counts
genusSlimmed <- plantPPE[complete.cases(plantPPE[,c("genus", "site", "age", "countPPE")]),
                                        c("genus", "site", "age", "countPPE")]
colnames(genusSlimmed)[colnames(genusSlimmed) == "countPPE"] = "count"

genus.novel <- neotoma.novelty(dataset = genusSlimmed,
                             ssmat.type = "abund",
                             bins = seq(-100,max(genusSlimmed$age), 200),
                             rich.cutoff = c(100, 5000),
                             age.limits = c(-150, 25100),
                             taxon.res = "genus",
                             bin.cutoff = 10,
                             taxa.cutoff = 5,
                             novel.alpha = 0.05,
                             novel.metric = "bray",
                             sqrt.mat=TRUE)

saveRDS(genus.novel, "./outputs/genus neotoma novelty.rds")

# ---------------------------- ####
# NOVEL TRAJECTORY CALCULATION ####

plant.record.df <- readRDS("./rawdata/processed_family_records.rds")
all.novel <- readRDS("./outputs/all neotoma novelty.rds")
site.df = read.csv("./rawdata/siteDf.csv")


fullRange <- c(3000,19000)

# NEED TO ADD RAREFACTION ETC TO RANDOM SUBSETS

novelTrajectories <- novel.trajectoryV4(novel.object = all.novel,
                                        ssmat.type = "prop.ssmats",
                                        dissim.method = "bray",
                                        novel.crop = fullRange)

saveRDS(novelTrajectories, date.wrap("./outputs/novelTrajectoryList",".rds"))

# ---------------------------- ####
# ####
# ------------------------ ####
# MODELLING PROCESS        ####
# ------------------------ ####

plant.record.df <- readRDS("./rawdata/processed_family_records.rds")
site.df = read.csv("./rawdata/siteDf.csv")
all.novel <- readRDS("./outputs/all neotoma novelty.rds")
novelTrajectories = readRDS("./outputs/novelTrajectoryList 2024-02-12.rds")

# IMPORT CLIMATE DATA ####
#           Convert temp records onto same scale ####

mod.env.data1 <- read.csv("./rawdata/climate/temp12k_allmethods_percentiles.csv")

cor.factor <- mean(mod.env.data1$global_median[mod.env.data1$age <= 11500 &
                                       mod.env.data1$age >= 6500])

old.env.data <- read.csv("./rawdata/climate/Shakun2012_retreat_temp.csv")
old.env.data$temp = old.env.data$temp + cor.factor
old.env.data$age = old.env.data$age * 1000

# Interpolate temp for each 200 year bin from reconstructions
mod.env.sub <- mod.env.data1[,c("ages","global_median")]
colnames(mod.env.sub) <- c("age", "temp")
comb.env.data <- rbind(old.env.data[old.env.data$age>max(mod.env.data1$ages),1:2], 
                       mod.env.sub)

temp.gam <- gam(temp ~ s(age, bs="cr", k=200), data=comb.env.data)
summary(temp.gam)

pred.df <- data.frame(age = seq(-100, 22000, 50))
temp.df <- cbind(pred.df,
                 as.data.frame(predict(temp.gam, newdata=pred.df, se.fit=TRUE)))

# modern comparison
temp.data <-  bin.env.data(env.data=temp.df,
                           bin.width=200,
                           lims=c(-100,23000),
                           env.var = "fit")
colnames(temp.data) <- gsub("env", "temp", colnames(temp.data))

write.csv(temp.data, "./outputs/interpolated temp data.csv")

#           local climate conditions ####

abs<-stack("./rawdata/climate/trace.01-36.22000BP.clm2.TSA.22000BP_decavg_400BCE.nc",
               varname="TSA")

extent(abs) = extent(c(0,360,-90,90))

# convert to celcius
abs = abs - 273.15

# dims are 1990AD (-50ybp) to 22000ybp in 10 year increments (first layer is oldest in raster)
traceBins = rev(seq(-40, 21990, 10))
neoBins = seq(-100,25000, 200)

targCoords = site.df[,c("long","lat")]
targCoords$long = ifelse(targCoords$long < 0, 360 + targCoords$long, targCoords$long)
colnames(targCoords) = c("x","y")

# get cell entry and try direct extract
dirEx = raster::extract(abs, targCoords)#[,-1] # remove ID

naCols = which(rowSums(is.na(dirEx))>0)
# custom function to extract adjacent temp cells for those on boundaries (e.g.,
#eastern American coast)

naVals = t(sapply(naCols, function(n){
  
  tempCoords = targCoords[n,]
  
  cellCoord = cellFromXY(abs, as.matrix(tempCoords))
  cellAdj = adjacent(abs, cells=cellCoord, directions="8")
  
  adjVals = raster::extract(abs, cellAdj[,2])
  adjVals = colMeans(adjVals, na.rm=TRUE)
  
  return(adjVals)
  
}))

traceMat = dirEx
traceMat[naCols,] = naVals

colnames(traceMat) = traceBins

# now iterate over each set of 200 and group by ID and take mean
traceCut = cut(traceBins, breaks=neoBins)
levels(traceCut) = neoBins[-length(neoBins)] + 0.5 * diff(neoBins[1:2])

traceMean = t(apply(traceMat, 1, function(x){tapply(x, traceCut, mean, na.rm=TRUE)}))
colnames(traceMean) = levels(traceCut)
traceMean = as.data.frame(traceMean)
traceMean$siteid = site.df$siteid

traceLong = as.data.frame(pivot_longer(traceMean, cols=colnames(traceMean)[colnames(traceMean) != "siteid"],
                                       names_to = "bin", values_to = "localTemp"))
summary(traceLong$localTemp)

# how many sites don't fall in raster cells?
traceSub = traceLong[traceLong$bin <= 22000,]
traceNa = unique(traceSub$site.id[is.na(traceSub$localTemp)])

# PREDICTIVE MODELS ####
#           Data prep ####
env.lags <- cbind(c(1,seq(5,25,5)),
                  c(1,seq(5,25,5)))

lagList <- lapply(c(1:nrow(env.lags)), function(n){
  
tempDf <- calcLag(novel.list = all.novel,
          env.data = temp.df,
          env.var = "fit",
          local.env.data = traceLong,
          global.lag=env.lags[n,1],
          local.lag=env.lags[n,2])

tempDf$binSite <- paste0(tempDf$site, ".", tempDf$bin)
return(tempDf)
})
lagDf <- cbind(lagList[[1]][,c(1,2,5)],
               do.call("cbind", lapply(lagList, function(x){
                 x[,3:4]
               })))

lagCor = cor(lagDf[,-(1:3)], use="complete.obs")
ifelse(abs(lagCor) > 0.3, round(lagCor, 3), 0)

post200 <- do.call("rbind", lapply(1:length(list(novelTrajectories)), function(n){
  x <- novelTrajectories$novel.traj.list[[3]]
  x<- x[x$time.since.novel <= 1000,]
  return(x)
}))

comb.df <- do.call("rbind", all.novel$novel)
comb.df <- comb.df[comb.df$cat == "novel",]
comb.df$novelID <- paste0(comb.df$site, ".", comb.df$bins)

post200 <- merge(post200, comb.df[,c("novelID", "cat.bef", "seq.dist", "raw.min.dist",
                                     "seq.exp", "min.exp", "min.p", "bin.lag")],
                 by.x="novelID", by.y="novelID",
                 all.x=TRUE, all.y=FALSE, sort=FALSE)

post200 <- merge(post200, temp.data,
                 by.x="bin", by.y="bin",
                 all.x=TRUE, all.y=FALSE, sort=FALSE)

traceLong$siteBin = paste0(traceLong$siteid, ".", traceLong$bin)

post200 <- merge(post200, traceLong[,-(1:2)],
                 by.x="novelID", by.y="siteBin",
                 all.x=TRUE, all.y=FALSE, sort=FALSE)

post200 <- merge(post200, lagDf[,-(1:2)],
                 by.x="novelID", by.y="binSite",
                 all.x=TRUE, all.y=FALSE, sort=FALSE)

post200 <- merge(post200, plant.record.df[!duplicated(plant.record.df$site),c("site", "lat", "long")],
                 by.x="site", by.y="site", all.x=TRUE, all.y=FALSE, sort=FALSE)

post200 <- post200[post200$novel.time <= 19000,]

# diversity change vars
post200$novelH1L <- post200$novelH1
post200$H1L <- post200$H1
post200$preNovH1L <- post200$preNovH1

post200$deltaH1 <- post200$novelH1L - post200$preNovH1L
post200$deltaH1abs <- abs(post200$deltaH1)

post200$deltaPost <- post200$H1L - post200$novelH1L
post200$deltaPostabs <- abs(post200$deltaPost)

post200$tsProp <- post200$bin.n / post200$tslength

post200$min.p <- 1-post200$min.p

post200$novAbundPropAbs <- ifelse(post200$novAbundProp > 1,
                                  post200$novAbundProp,
                                  1/post200$novAbundProp)
plot(post200$novAbundPropAbs ~ post200$novAbundProp)

tempRot <- rotate.data(post200$dP, post200$dN, -45)

post200$dMag <- tempRot[,3]
post200$dRatio <- tempRot[,4]

scaleVars = c("novelH1L", "deltaH1",
              "deltaPost", "deltaH1abs", "deltaPostabs",
              "min.exp",
              "bin.lag", "time.since.novel",
              "gamma", "tslength", "tsProp",
              "temp", "localTemp","min.p",
              colnames(post200)[grepl("Lag", colnames(post200))],
              "novAbundPropAbs",
              "novAbundProp")

# make sure we don't have any NA rows so we can conduct model comparisons
post200 <- post200[complete.cases(post200[,c("dMag","dRatio", scaleVars)]),]

scaleVars <- lapply(post200[, c("novelH1L", "deltaH1",
                                "deltaPost", "deltaH1abs", "deltaPostabs",
                                "min.exp",
                                "bin.lag", "time.since.novel",
                                "gamma", "tslength", "tsProp",
                                "temp", "localTemp","min.p",
                                colnames(post200)[grepl("Lag", colnames(post200))],
                                "novAbundPropAbs",
                                "novAbundProp")], function(x){
                                  x <- scale(x)
                                  return(list(c(attr(x, "scaled:center"),
                                                attr(x, "scaled:scale")),
                                              as.vector(x)))
                                })

scaleDf <- t(sapply(scaleVars, function(x){x[[1]]}))
scaleVars <- as.data.frame(do.call("cbind", lapply(scaleVars, function(x){x[[2]]})))
colnames(scaleVars) <- paste0(colnames(scaleVars), "S")

post200 <- cbind(post200, scaleVars)

#           forward model selection ####

post200$dMagRaw <- post200$dMag
post200$dMag <- log(post200$dMag)
post200$novelID = as.factor(post200$novelID)

post200$gaProp = post200$novelH0 / post200$gamma 
post200$gaPropS = scale(post200$gaProp)

write.csv(post200, "./outputs/post200.csv")

# intercept only
dMint <- lmer(dMag ~ 1 + (1|novelID), data=post200)
dRint <- lmer(dRatio ~ 1 + (1|novelID), data=post200)

dMintSite = lmer(dMag ~ 1 + (1|site/novelID), data=post200)
dRintSite = lmer(dRatio ~ 1 + (1|site/novelID), data=post200)

compare_performance(dMint, dMintSite)
compare_performance(dRint, dRintSite)

# novelty covariates
dMcov <- update(dMint, .~. + bin.lagS + tsPropS + tslengthS + gammaS + gaPropS + novAbundPropS + novAbundPropAbsS)
dRcov <- update(dRint, .~. + bin.lagS + tsPropS + tslengthS + gammaS + gaPropS + novAbundPropS + novAbundPropAbsS)

# distance decay
dMdecay <- update(dMcov, .~. + time.since.novelS)
dRdecay <- update(dRcov, .~. + time.since.novelS)

# novelty expectations
dMnovelty <- update(dMdecay, .~. + min.pS)
dRnovelty <- update(dRdecay, .~. + min.pS)
  
# climate

# iterate through lag options
lagGrid <- expand.grid(colnames(post200)[grepl("globalLag", colnames(post200)) & grepl("S", colnames(post200))],
                       colnames(post200)[grepl("localLag", colnames(post200)) & grepl("S", colnames(post200))])
climAIC <- do.call("rbind", lapply(1:nrow(lagGrid), function(n){
  print(n)
  post200$tempGlobal <- post200[,as.character(lagGrid[n,1])]
  post200$tempLocal <- post200[,as.character(lagGrid[n,2])]
  
  post200$tempGlobalA = scale(abs(post200[,gsub("S","",as.character(lagGrid[n,1]))]))
  post200$tempLocalA = scale(abs(post200[,gsub("S","",as.character(lagGrid[n,2]))]))
  
  return(data.frame(globalLag = lagGrid[n,1],
                    localLag = lagGrid[n,2],
                    dM = AIC(update(dMnovelty, .~. + tempGlobal + tempLocal)),
                    dR = AIC(update(dRnovelty, .~. + tempGlobal + tempLocal)),
                    dMa = AIC(update(dMnovelty, .~. + tempGlobalA + tempLocalA)),
                    dRa = AIC(update(dRnovelty, .~. + tempGlobalA + tempLocalA))))
  
}))

plot(climAIC[,3] ~ climAIC[,4], type="n")
text(climAIC[,3] ~ climAIC[,4])
# best fit is for both M and R models (actually like 0.25AIC worse for dR...)
climVars = climAIC[which.min(rowSums(climAIC[,3:4])),1:2]

dMclim = update(dMnovelty, as.formula(paste0(".~. + ", climVars[1,1], " + ", climVars[1,2])))
dRclim = update(dRnovelty, as.formula(paste0(".~. + ", climVars[1,1], " + ", climVars[1,2])))

# novel diversity DIRECTIONAL
dMdivNoPost <- update(dMclim, .~. + (novelH1LS + deltaH1S)^2)
dRdivNoPost <- update(dRclim, .~. + (novelH1LS + deltaH1S)^2)

# novel diversity ABSOLUTE
dMdivNoPostabs <- update(dMclim, .~. + (novelH1LS + deltaH1absS)^2)
dRdivNoPostabs <- update(dRclim, .~. + (novelH1LS + deltaH1absS)^2)

# post-novel DIRECTIONAL (novel DIRECTIONAL)
dMdivDirDir <- update(dMclim, .~. + (novelH1LS + deltaH1S + deltaPostS)^2)
dRdivDirDir <- update(dRclim, .~. + (novelH1LS + deltaH1S + deltaPostS)^2)

# post-novel ABSOLUTE (novel DIRECTIONAL)
dMdivDirAbs <- update(dMclim, .~. + (novelH1LS + deltaH1S + deltaPostabsS)^2)
dRdivDirAbs <- update(dRclim, .~. + (novelH1LS + deltaH1S + deltaPostabsS)^2)

# post-novel DIRECTIONAL (novel ABSOLUTE)
dMdivAbsDir <- update(dMclim, .~. + (novelH1LS + deltaH1absS + deltaPostS)^2)
dRdivAbsDir <- update(dRclim, .~. + (novelH1LS + deltaH1absS + deltaPostS)^2)

# post-novel ABSOLUTE (novel ABSOLUTE)
dMdivAbsAbs <- update(dMclim, .~. + (novelH1LS + deltaH1absS + deltaPostabsS)^2)
dRdivAbsAbs <- update(dRclim, .~. + (novelH1LS + deltaH1absS + deltaPostabsS)^2)

# div climate interactions?
dRdivClimInt <- update(dRdivDirDir, as.formula(paste0(".~. + (novelH1LS + deltaH1S + deltaPostS) * (", climVars[1,1], " + ", climVars[1,2], ")")))
dMdivClimInt <- update(dMdivAbsAbs, as.formula(paste0(".~. + (novelH1LS + deltaH1absS + deltaPostabsS) * (", climVars[1,1], " + ", climVars[1,2], ")")))

dMcomp <- compare_performance(dMint, dMcov, dMdecay, dMnovelty, dMclim,
                              dMdivNoPost, dMdivNoPostabs,
                              dMdivDirDir, dMdivAbsDir,dMdivDirAbs,dMdivAbsAbs)

write.csv(dMcomp, "./outputs/dMcomp.csv")

dRcomp <- compare_performance(dRint, dRcov, dRdecay, dRnovelty, dRclim, 
                              dRdivNoPost, dRdivNoPostabs,
                              dRdivDirDir, dRdivAbsDir,dRdivDirAbs,dRdivAbsAbs)

write.csv(dRcomp, "./outputs/dRcomp.csv")

write.csv(summary(dMdivAbsAbs)$coefficients, "./outputs/dMdivAbsAbs.csv")
write.csv(summary(dRdivDirDir)$coefficients, "./outputs/dRdivDirDir.csv")

#           model assumption tests ####

turnRes <- simulateResiduals(dMdivAbsAbs)
plot(turnRes)

turnTests <- modelDiagTests(dMdivAbsAbs,
                            time=post200$bin,
                            coords=post200[,c("lat","long")],
                            spat.iter=999)
turnTests$unif
turnTests$disp
summary(turnTests$tauto$statistic)
summary(turnTests$tauto$p <= 0.05)
summary(turnTests$sauto$observed)
summary(turnTests$sauto$p <= 0.05)

persTests <- modelDiagTests(dRdivDirDir,
                            time=post200$bin,
                            coords=post200[,c("lat","long")],
                            spat.iter=999)
persTests$unif
persTests$disp
summary(persTests$tauto$statistic)
summary(persTests$tauto$p <= 0.05)
summary(persTests$sauto$observed)
summary(persTests$sauto$p <= 0.05)

#           random non-novel models ####

rNgamList <- lapply(1:999, function(n){

  print(n)
  
  randX <- do.call("rbind", lapply(1, function(n){
    
    x <- novelTrajectories
    rand200 <- x$rand.traj.list[x$rand.traj.list$novel.time <= 19000,]
    # pick random non-novel points equal to novel sampling
    novLength <- length(unique(post200$novelID))
    randIDs <- sample(unique(rand200$novelID), novLength)
    
    return(rand200[rand200$novelID %in% randIDs & rand200$time.since.novel <= 1000,])
    
  }))
  print(nrow(randX))
  
  tempRot <- rotate.data(randX$dP, randX$dN, -45)
  
  randX$dMag <- tempRot[,3]
  randX$dRatio <- tempRot[,4]
  randX$dMag <- log(randX$dMag)
  
  comb.df <- do.call("rbind", all.novel$novel)
  comb.df$novelID <- paste0(comb.df$site, ".", comb.df$bins)
  
  randX <- merge(randX, comb.df[,c("novelID", "cat.bef", "seq.dist", "raw.min.dist",
                                   "seq.exp", "min.exp", "min.p", "bin.lag")],
                 by.x="novelID", by.y="novelID",
                 all.x=TRUE, all.y=FALSE, sort=FALSE)
  
  randX <- merge(randX, lagDf[,-(1:2)],
                   by.x="novelID", by.y="binSite",
                   all.x=TRUE, all.y=FALSE, sort=FALSE)
  
  randX$tsProp <- randX$bin.n / randX$tslength
  
  randX$novelH1L <- log(randX$novelH1+1)
  randX$H1L <- log(randX$H1+1)
  randX$preNovH1L <- log(randX$preNovH1+1)
  
  randX$deltaH1 <- randX$novelH1L - randX$preNovH1L
  randX$deltaH1abs <- abs(randX$deltaH1)
  
  randX$deltaPost <- randX$H1L - randX$novelH1L
  randX$deltaPostabs <- abs(randX$deltaPost)
  
  randX$gaProp = randX$novelH0 / randX$gamma 
  randX$gaPropS = scale(randX$gaProp)
  
  randX$min.p <- 1-randX$min.p
  
  randX$novAbundPropAbs <- ifelse(randX$novAbundProp > 1,
                                  randX$novAbundProp,
                                    1/randX$novAbundProp)
  
  randscaleVars <- lapply(randX[, c("novelH1L", "deltaH1", "deltaH1abs", "deltaPost", "deltaPostabs",
                                    "min.exp", "novAbundProp", "novAbundPropAbs",
                                    "min.p",
                                    "bin.lag", "time.since.novel",
                                    "gamma", "tslength", "tsProp",
                                    "localLag15", "globalLag20")], function(x){
                                  x <- scale(x)
                                  return(list(c(attr(x, "scaled:center"),
                                                attr(x, "scaled:scale")),
                                              as.vector(x)))
                                })
  
  scaleVars <- as.data.frame(do.call("cbind", lapply(randscaleVars, function(x){x[[2]]})))
  colnames(scaleVars) <- paste0(colnames(scaleVars), "S")
  
  randX <- cbind(randX, scaleVars)
  
  randX$timeCat <- as.factor(randX$time.since.novel)
  randX$site <- as.factor(randX$site)

  randX$novelID <- as.factor(randX$novelID)
  
  randX <- randX[complete.cases(scaleVars),]
  
  # model selection (for performance)
  # intercept only
  
  # intercept only
  rMint <- lmer(dMag ~ 1 + (1|novelID), data=randX)
  rRint <- lmer(dRatio ~ 1 + (1|novelID), data=randX)
  
  # novelty covariates
  rMcov <- update(rMint, .~. + bin.lagS + tsPropS + tslengthS + gammaS + gaPropS + novAbundPropS + novAbundPropAbsS)
  rRcov <- update(rRint, .~. + bin.lagS + tsPropS + tslengthS + gammaS + gaPropS + novAbundPropS + novAbundPropAbsS)
  
  # distance decay
  rMdecay <- update(rMcov, .~. + time.since.novelS)
  rRdecay <- update(rRcov, .~. + time.since.novelS)
  
  # novelty expectations
  rMnovelty <- update(rMdecay, .~. + min.pS)
  rRnovelty <- update(rRdecay, .~. + min.pS)
  
  # climate
  rMclim = update(rMnovelty, as.formula(paste0(".~. + ", climVars[1,1], " + ", climVars[1,2])))
  rRclim = update(rRnovelty, as.formula(paste0(".~. + ", climVars[1,1], " + ", climVars[1,2])))
  
  # novel diversity DIRECTIONAL
  rMdivNoPost <- update(rMclim, .~. + (novelH1LS + deltaH1S)^2)
  rRdivNoPost <- update(rRclim, .~. + (novelH1LS + deltaH1S)^2)
  
  # novel diversity ABSOLUTE
  rMdivNoPostabs <- update(rMclim, .~. + (novelH1LS + deltaH1absS)^2)
  rRdivNoPostabs <- update(rRclim, .~. + (novelH1LS + deltaH1absS)^2)
  
  # post-novel DIRECTIONAL (novel DIRECTIONAL)
  rMdivDirDir <- update(rMclim, .~. + (novelH1LS + deltaH1S + deltaPostS)^2)
  rRdivDirDir <- update(rRclim, .~. + (novelH1LS + deltaH1S + deltaPostS)^2)
  
  # post-novel ABSOLUTE (novel DIRECTIONAL)
  rMdivDirAbs <- update(rMclim, .~. + (novelH1LS + deltaH1S + deltaPostabsS)^2)
  rRdivDirAbs <- update(rRclim, .~. + (novelH1LS + deltaH1S + deltaPostabsS)^2)
  
  # post-novel DIRECTIONAL (novel ABSOLUTE)
  rMdivAbsDir <- update(rMclim, .~. + (novelH1LS + deltaH1absS + deltaPostS)^2)
  rRdivAbsDir <- update(rRclim, .~. + (novelH1LS + deltaH1absS + deltaPostS)^2)
  
  # post-novel ABSOLUTE (novel ABSOLUTE)
  rMdivAbsAbs <- update(rMclim, .~. + (novelH1LS + deltaH1absS + deltaPostabsS)^2)
  rRdivAbsAbs <- update(rRclim, .~. + (novelH1LS + deltaH1absS + deltaPostabsS)^2)

  # Final models (for coefficients)  
  randdMgam <- rMdivAbsAbs
  
  randdMcoef <- summary(randdMgam)$coefficients
  
  randdRgam <- rRdivDirDir
  
  randdRcoef <- summary(randdRgam)$coefficients
  
 # Diversity predictions
  pred.dfM <- expand.grid(deltaPostabsS = c(quantile(post200$deltaPostabsS, c(0.1)),
                                            0,
                                            quantile(post200$deltaPostabsS, c(0.9))),
                          deltaH1absS = c(quantile(post200$deltaH1absS, c(0.1)),
                                          0,
                                          quantile(post200$deltaH1absS, c(0.9))),
                          novelH1LS = rev(quantile(post200$novelH1LS, c(0.1, 0.5, 0.90))))
  
  predMM <- matrix(0, nrow=nrow(pred.dfM), ncol=ncol(randdMgam@frame), dimnames=list(NULL, colnames(randdMgam@frame)))
  predMM <- predMM[,!colnames(predMM) %in% c("dMag", "novelID")]
  predMM <- as.data.frame(predMM)
  predMM[,colnames(pred.dfM)] = pred.dfM
  
  divMpreds <- predict(randdMgam, newdata=predMM, re.form=NA)
  divMpreds <- exp(divMpreds)
  
  pred.df <- expand.grid(deltaPostS = c(quantile(randX$deltaPostS, c(0.1)),
                                        (0 - scaleDf["deltaPost",1]) / scaleDf["deltaPost",2],
                                        quantile(randX$deltaPostS, c(0.9))),
                         deltaH1S = c(quantile(randX$deltaH1S, c(0.1)),
                                      (0 - scaleDf["deltaH1",1]) / scaleDf["deltaH1",2],
                                      quantile(randX$deltaH1S, c(0.9))),
                         novelH1LS = rev(quantile(randX$novelH1LS, c(0.1, 0.5, 0.90))))
  
  predMM <- matrix(0, nrow=nrow(pred.df), ncol=ncol(randdRgam@frame), dimnames=list(NULL, colnames(randdRgam@frame)))
  predMM <- predMM[,!colnames(predMM) %in% c("dRatio", "novelID")]
  predMM <- as.data.frame(predMM)
  predMM[,colnames(pred.df)] = pred.df
  
  divRpreds <- predict(randdRgam, newdata=predMM, re.form=NA)
  
  return(list(Mcoef = randdMcoef,
              Rcoef = randdRcoef,
              Mperf = compare_performance(rMint, rMcov, rMdecay, rMnovelty, rMclim, rMdivNoPostabs, randdMgam),
              Rperf = compare_performance(rRint, rRcov, rRdecay, rRnovelty, rRclim, rRdivNoPost, randdRgam),
              Mpreds = divMpreds,
              Rpreds = divRpreds,
              Mfull = randdMgam,
              Rfull = randdRgam))
              
})

library(abind)

dMcoef <- summary(dMdivAbsAbs)$coefficients
dRcoef <- summary(dRdivDirDir)$coefficients

randMcoef <- abind(lapply(rNgamList, function(x){x$Mcoef}), along=3, force.array=TRUE)
randRcoef <- abind(lapply(rNgamList, function(x){x$Rcoef}), along=3, force.array=TRUE)

randMr2 <- do.call("rbind", lapply(rNgamList, function(x){x$Mperf$R2_marginal}))
randRr2 <- do.call("rbind", lapply(rNgamList, function(x){x$Rperf$R2_marginal}))

fullMr2 = do.call("rbind", lapply(rNgamList, function(x){performance(x$Mfull)}))
summary(fullMr2$R2_marginal)
summary(fullMr2$R2_conditional)

fullRr2 = do.call("rbind", lapply(rNgamList, function(x){performance(x$Rfull)}))
summary(fullRr2$R2_marginal)
summary(fullRr2$R2_conditional)

# ----- ####
# PLOTS ####
# ----- ####
#           color ramps ####

turnRamp <- colorRamp(c("grey90", "grey50"))
persRamp <- colorRamp(c("purple", "grey80", "orange"))

#           Raw observations ####

# null distribution of Turn and Pers
rTPlist <- lapply(1:length(rNgamList), function(n){
  
  print(n)
  
  randX <- do.call("rbind", lapply(1, function(n){
    
    x <- novelTrajectories
    rand200 <- x$rand.traj.list[x$rand.traj.list$novel.time <= 19000,]
    # pick random non-novel points equal to novel sampling
    novLength <- length(unique(post200$novelID))
    randIDs <- sample(unique(rand200$novelID), novLength)
    
    return(rand200[rand200$novelID %in% randIDs & rand200$time.since.novel <= 1000,])
    
  }))
  print(nrow(randX))
  
  tempRot <- rotate.data(randX$dP, randX$dN, -45)
  
  randX$dMag <- tempRot[,3]
  randX$dRatio <- tempRot[,4]
  randX$dMag <- log(randX$dMag)
  
  return(cbind(randX$dMag, randX$dRatio))
})

# gam of mag to ratio relationship
library(gamm4)
post200$time.since.novelFact = as.factor(post200$time.since.novel)

MRgam = lmer(dRatio ~ dMag * time.since.novelFact + (1|novelID), data=post200)
ratPred = data.frame(dMag=rep(seq(min(post200$dMag), max(post200$dMag), len=200), 5),
                     time.since.novelFact=rep(levels(post200$time.since.novelFact),
                                              each=200))

pred.df = cbind(ratPred, mer.ci(MRgam, ratPred, sims=99))

pdf(date.wrap("./plots/rawRotation", ".pdf"), width=7, height=4, useDingbats = FALSE)

split.screen(rbind(c(0.1,0.6075,0.1,0.99),
                   c(0.675,0.99,0.61,0.99),
                   c(0.675,0.99,0.1,0.48)))

screen(1)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1, mgp=c(3,0.5,0))

plot(x=NULL, y=NULL, ylim=quantile(post200$dRatio, probs=c(0.01,0.99)), 
     xlim=c(0,1.1), xlab="", ylab="",
     xaxt="n", yaxs="i", xaxs="i")
axis(side=1, mgp=c(3,0.1,0))
mtext(side=1, line=1, text="Turnover")
mtext(side=2, line=1.75, text="Persistence", las=0)
abline(h=0, col="grey", lty="31")

# color res
colRes <- 0.001
xSeq <- seq(0, max(post200$dMag)+colRes, len=2000)
ySeq <- seq(min(post200$dRatio), max(post200$dRatio), len=2000)
topGrid <- colorRampPalette(c("white", "orange"))(length(ySeq))
botGrid <- colorRampPalette(c("purple", "black"))(length(xSeq))
colGrid <- mapply(x=topGrid, y=botGrid, function(x,y){
  colorRampPalette(c(x,y))(length(xSeq))
})
gridX <- match(round(post200$dMag,3), round(xSeq,3))
gridY <- match(round(post200$dRatio,3), round(ySeq,3))

points(post200$dRatio ~ post200$dMagRaw, pch=16, cex=0.3,
       col="grey85")#       col=colGrid[cbind(gridX,gridY)], lwd=0.5)

rampRange = max(abs(range(pred.df$fit)))

pred.df$predRamp = unitScale(pred.df$fit, custMin= -rampRange, custMax= rampRange)

predSplit = split(pred.df, f=pred.df$time.since.novelFact)
lapply(predSplit, function(x){

  obsSum = range(post200$dMagRaw[post200$time.since.novel == x$time.since.novelFact[1]])
  x$dMag = exp(x$dMag)
  x = x[x$dMag >= obsSum[1] & x$dMag <= obsSum[2],]
  
sapply(2:nrow(x), function(n){
  print(n)
  rampCol = colorRamp(c("purple","grey70","orange"))(x$predRamp[n])/255
  
  lines(x$fit[(n-1):n] ~ x$dMag[(n-1):n], lwd=2,
        col=rgb(rampCol[1], rampCol[2], rampCol[3]))
  
  polygon(x=c(x$dMag[(n-1):n], rev(x$dMag[(n-1):n])),
          y=c(x$upper[(n-1):n], rev(x$lower[(n-1):n])),
          col=rgb(rampCol[1], rampCol[2], rampCol[3], 0.5),
          border=NA)
  
})

  text(x=rev(x$dMag)[1],
       y=rev(x$fit)[1],
       labels=x$time.since.novelFact[1], pos=4)

})

text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.05, "y"),
     labels=expression("R"^2*" = "), adj=0)

text(x=relative.axis.point(0.13, "x"),
     y=relative.axis.point(0.04, "y"),
     labels=sprintf("%.3f", performance(MRgam)$R2_marginal), adj=0)

text(x=relative.axis.point(0.015, "x"),
     y=relative.axis.point(0.975, "y"),
     labels="(A)", font=2, adj=0)

box()
close.screen(1)

screen(2)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1, mgp=c(3,0.5,0))

turnDens <- lapply(sort(unique(post200$time.since.novel)), function(x){
  subPost = post200[post200$time.since.novel == x,]                 
                density(subPost$dMagRaw,
                    from=min(post200$dMagRaw),
                    to=max(post200$dMagRaw))
})

plot(NULL, xlim=c(0,1), ylim=c(0,5), xlab="", ylab="", xaxt="n")

sapply(1:length(turnDens), function(xN){
  x = turnDens[[xN]]
  unitX = unitScale(x$x, custMin = min(post200$dMagRaw), custMax = max(post200$dMagRaw))
  sapply(2:length(x$x), function(n){
  lines(x=x$x[(n-1):n], y=x$y[(n-1):n], col=rgb(turnRamp(unitX[n])/255), lwd=2)
  })
  text(y=max(x$y), x=x$x[which.max(x$y)], pos=4,
       labels=seq(200,1000,200)[xN], col="grey50")
})

axis(side=1, mgp=c(3,0.1,0), at=c(0,0.25,0.5,0.75,1,1.25))

rTT = t(sapply(rTPlist, function(x){density(exp(x[,1]), from=0, to=sqrt(2))$y}))
rTTm <- apply(rTT, 2, mean)
lines(rTTm ~ seq(0, sqrt(2), len=length(rTTm)), lwd=1.5, lty="31")

mtext(side=1, text="Turnover", line=1)
mtext(side=2, text="Density", las=0, line=1.5)

text(x=c(0.25, 0.45), y=c(4.5, 2), offset=0.25,
     labels=c("Null", "Observed"), pos=4)

text(x=relative.axis.point(0.025, "x"),
     y=relative.axis.point(0.925, "y"),
     labels="(B)", font=2, adj=0)
close.screen(2)

screen(3)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1, mgp=c(3,0.5,0))

persDens <- lapply(sort(unique(post200$time.since.novel)), function(x){
  subPost = post200[post200$time.since.novel == x,]                 
  density(subPost$dRatio,
          from=min(post200$dRatio),
          to=max(post200$dRatio))
})

plot(NULL, xlim=c(-0.5,0.5), ylim=c(0,12), xlab="", ylab="", xaxt="n")

sapply(1:length(persDens), function(xN){
  x = persDens[[xN]]
  unitX = unitScale(x$x, custMin = min(post200$dRatio), custMax = max(post200$dRatio))
  sapply(2:length(x$x), function(n){
    lines(x=x$x[(n-1):n], y=x$y[(n-1):n], col=rgb(persRamp(unitX[n])/255), lwd=2)
  })
  text(y=max(x$y), x=x$x[which.max(x$y)], pos=4,
       labels=seq(200,1000,200)[xN], col="grey50")
  return(unitX)
})

axis(side=1, mgp=c(3,0.1,0))

rTP = t(sapply(rTPlist, function(x){density(x[,2], from=-sqrt(0.5), to=sqrt(0.5))$y}))
rTPm <- apply(rTP, 2, mean)
lines(rTPm ~ seq(-sqrt(0.5), sqrt(0.5), len=length(rTPm)), lwd=1.25, lty="31")

mtext(side=1, text="Persistence", line=1)
mtext(side=2, text="Density", las=0, line=1.5)

text(x=c(0.025, 0.15), y=c(10, 3), labels=c("Null", "Observed"), pos=4, offset=0.25)

text(x=relative.axis.point(0.025, "x"),
     y=relative.axis.point(0.925, "y"),
     labels="(C)", font=2, adj=0)
close.screen(3)

dev.off()

#           model performance ####

Mrank <- c(1,2,3,4,5,rep(6,2), rep(7,4))
dMcomp <- do.call("rbind", lapply(split(dMcomp, f=Mrank), function(x){
  x[which.min(x$AIC),]
}))
dMdAIC = dMcomp$R2_marginal

dRcomp <- do.call("rbind", lapply(split(dRcomp, f=Mrank), function(x){
  x[which.min(x$AIC),]
}))
dRdAIC = dRcomp$R2_marginal

pdf(date.wrap("./plots/newModels", ".pdf"), height=4.5, width=6.5, useDingbats=FALSE)

par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)

split.screen(rbind(c(0.3,0.625,0.55,0.98),
                   c(0.625,0.95,0.55,0.98),
                   c(0.3,0.625,0.15,0.55),
                   c(0.625,0.95,0.15,0.55)))

# Magnitude #
screen(1)

plot(x=NULL, y=NULL, 
     xlim=c(0,0.35), 
     ylim=rev(c(0,6.25)), xaxs="i",
     xlab="", ylab="", yaxt="n", xaxt="n")
axis(side=1, mgp=c(3,0.2,0), at=seq(0,0.3,0.1), labels=NA)
#segments(x0=0, x1=0, y0=par("usr")[3], y1=relative.axis.point(0.9,"y"), lty="31")

randMr2d = t(apply(randMr2, 1, diff))

rect(xleft=apply(randMr2, 2, quantile, 0.025),
     xright=apply(randMr2, 2, quantile, 0.975),
     ytop=0:6 + 0.25,
     ybottom = 0:6 - 0.25,
     col="grey70")
points(y=1:6, x=dMdAIC[-1], pch=16)

axis(side=2, at=1:6,
     labels=paste0("(",2:7,")"), las=1)
par(lheight=0.7)
mtext(side=2, at=relative.axis.point(0.95, "y"),
      text="Model\nstage", line=0.25, las=1)

par(xpd=NA)
axis(side=2, at=c(0.75,6.25), line=2, tcl=0.25, labels=NA)
par(xpd=FALSE)

axis(side=2, at=1:6, line=2,
     labels=c("Covariates", "Time decay",
              "Novelty extent", "Climate change",
              "Novel\ndiversity change",
              "Post-novel\ndiversity change"), las=1)
par(lheight=1)

text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.95, "y"),
     labels=expression(bold("(A)")*" Turnover"),
     font=2, adj=0)

close.screen(1)

screen(2)

plot(x=NULL, y=NULL, 
     xlim=c(0,0.18), 
     ylim=rev(c(0,6.25)), xaxs="i",
     xlab="", ylab="", yaxt="n", xaxt="n")
axis(side=1, mgp=c(3,0.2,0), at=seq(0,0.16,0.04), labels=NA)
segments(x0=0, x1=0, y0=par("usr")[3], y1=relative.axis.point(0.9,"y"), lty="31")

randMr2d = t(apply(randMr2, 1, diff))

rect(xleft=apply(randMr2d, 2, quantile, 0.025),
     xright=apply(randMr2d, 2, quantile, 0.975),
     ytop=1:6 + 0.25,
     ybottom = 1:6 - 0.25,
     col="grey70")
points(y=1:6, x=diff(dMdAIC), pch=16)

axis(side=2, at=1:6,
     labels=NA, las=1)

text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.95, "y"),
     labels=expression(bold("(B)")*" Turnover"),
     font=2, adj=0)

close.screen(2)

# Persistence #
screen(3)
plot(x=NULL, y=NULL, 
     xlim=c(0,0.35), 
     ylim=rev(c(0,6.25)), xaxs="i",
     xlab="", ylab="", yaxt="n", xaxt="n")
axis(side=1, mgp=c(3,0.2,0))
mtext(side=1, line=1, text="Model performance")
mtext(side=1, line=1.9, text=expression("(Marginal R"^2*")"))
                      
#segments(x0=0, x1=0, y0=par("usr")[3], y1=relative.axis.point(0.9,"y"), lty="31")

randRr2d = t(apply(randRr2, 1, diff))

rect(xleft=apply(randRr2, 2, quantile, 0.025),
     xright=apply(randRr2, 2, quantile, 0.975),
     ytop=0:6 + 0.25,
     ybottom = 0:6 - 0.25,
     col=colorRampPalette(c("purple","orange"))(3)[2])
points(y=1:6, x=dRdAIC[-1], pch=16)

axis(side=2, at=1:6,
     labels=paste0("(",2:7,")"), las=1)
par(lheight=0.7)
par(xpd=NA)
axis(side=2, at=c(0.75,6.25), line=2, tcl=0.25, labels=NA)
par(xpd=FALSE)

axis(side=2, at=1:6, line=2,
     labels=c("Covariates", "Time decay",
              "Novelty extent", "Climate change",
              "Novel\ndiversity change",
              "Post-novel\ndiversity change"), las=1)
par(lheight=1)


text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.925, "y"),
     labels=expression(bold("(C)")*" Persistence"),
     font=2, adj=0)

close.screen(3)

# Persistence #
screen(4)
plot(x=NULL, y=NULL, 
     xlim=c(0,0.18), 
     ylim=rev(c(0,6.25)), xaxs="i",
     xlab="", ylab="", yaxt="n", xaxt="n")
axis(side=1, mgp=c(3,0.2,0))

mtext(side=1, line=1, text="Model improvement")
mtext(side=1, line=1.9, text=expression("("*Delta*" Marginal R"^2*")"))

#segments(x0=0, x1=0, y0=par("usr")[3], y1=relative.axis.point(0.9,"y"), lty="31")

randRr2d = t(apply(randRr2, 1, diff))

rect(xleft=apply(randRr2d, 2, quantile, 0.025),
     xright=apply(randRr2d, 2, quantile, 0.975),
     ytop=1:6 + 0.25,
     ybottom = 1:6 - 0.25,
     col=colorRampPalette(c("purple","orange"))(3)[2])
points(y=1:6, x=diff(dRdAIC), pch=16)
axis(side=2, at=1:6, labels=NA, las=1)
text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.925, "y"),
     labels=expression(bold("(D)")*" Persistence"),
     font=2, adj=0)

legend(x=relative.axis.point(1, "x"), y=relative.axis.point(0.6, "y"),
       legend=c("Obs", "Null"),
       fill=c(NULL, "white"),
       border=c("white", "black"), bty="n",
       pch=c(16, NA), xjust=1, y.intersp=0.75)
rect(xleft=relative.axis.point(0.655, "x"),
     xright=relative.axis.point(1, "x"),
     ybottom=relative.axis.point(0.355, "y"),
     ytop=relative.axis.point(0.575, "y"))

close.screen(4)

close.screen(all.screens=TRUE)
dev.off()

#           coef plot ####

plot.order <- rownames(dMcoef)[-1]
plotLabels = plot.order 

plot.order <- c("bin.lagS",
                "tsPropS",
                "tslengthS",
                "gammaS",
                "gaPropS",
                "novAbundPropS",
                "novAbundPropAbsS",
                "time.since.novelS",
                "min.pS",
                "globalLag20S", "localLag15S",
                "novelH1LS", "deltaH1S", "novelH1LS:deltaH1S",
                "deltaPostS", "novelH1LS:deltaPostS", "deltaH1S:deltaPostS")

plotLabels <- c("Time lag",
                "Time series position",
                "Time series length",
                "Gamma diversity",
                "Gamma proportion",
                "Pollen variation",
                "Abs pollen variation",
                "Time since novelty",
                "Novelty magnitude",
                expression(Delta*"Global temp"),
                expression(Delta*"Regional temp"),
                "Novel state H1 (N)",
                expression("Novel "*Delta*"H1 ("*Delta*"N)"),
                expression(underline("N:"*Delta*"N")),
                expression("Post-novel "*Delta*"H1 ("*Delta*"pN)"),
                expression(underline("N':"*Delta*"N")),
                expression(underline(Delta*"N:"*Delta*"pN")))

pdf(date.wrap("./plots/predCoefPlot", ".pdf"), height=5.5, width=8.5, useDingbats=FALSE,
    bg="white")

xlims <- c(-0.2,0.4)
offset <- 0

split.screen(rbind(c(0.25,0.62,0.09,0.99),
                   c(0.62,0.99,0.09,0.99)))

screen(1)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0))

plot(x=NULL, y=NULL, xlim=c(-0.12,0.135), ylim=rev(c(0.25, length(plot.order))), yaxt="n", xlab="",
     ylab="", xaxt="n")

# dM rand coefs
dMrands <- randMcoef[,1,]
obsCoef <- dMcoef

rownames(obsCoef) = gsub('abs', "", rownames(obsCoef))
rownames(dMrands) = gsub('abs', "", rownames(dMrands))

# rearrange to right order
obsCoef <- obsCoef[match(plot.order, rownames(obsCoef)),]
dMrands <- dMrands[match(plot.order, rownames(dMrands)),]

boxWidth = 0.5

# background dividers
bDivs <- c(0,7,8,9,11,14,par("usr")[3])

rect(xleft=par("usr")[1], xright=par("usr")[2],
     ybottom=bDivs[-length(bDivs)]+0.5, ytop=bDivs[-1]+0.5,
     col=c("grey90","white"), border=NA)
text(y=bDivs[-length(bDivs)]+0.9,
     x=c(rep(relative.axis.point(0.96, "x"),3),
         rep(relative.axis.point(0.96, "x"),2)),
     labels=rev(c("(6) Post-novel diversity",
              "(5) Novel diversity",
              "(4) Climate",
              "(3) Novelty magnitude",
              "(2) Time decay",
              "(1) Covariates")),
     adj=1, font=2, col="grey70")
# text(y=bDivs[3]+0.75,
#      x=relative.axis.point(0.95, "x"),
#      labels=c("Novelty"),
#      adj=0, font=2, col="grey70")

segments(x0=0, x1=0, y0=par("usr")[3], y1=relative.axis.point(0.95,"y"),
         lty="31")
axis(side=1, mgp=c(3,0.1,0))

greyMean <- apply(dMrands, 1, mean)
greyMean <- unitScale(greyMean)

sapply(1:nrow(dMrands), function(n){
  print(n)
rect(xleft=quantile(dMrands[n,], prob=0.025),
     xright=quantile(dMrands[n,], prob=0.975),
     ybottom=n - 0.5 * boxWidth - offset,
     ytop = n  + 0.5 * boxWidth - offset,
     col=rgb(turnRamp(greyMean[n])/255))
})
  
segments(y0=1:nrow(obsCoef) - offset, y1=1:nrow(obsCoef) - offset,
         x0=obsCoef[,1] + 1.96 * obsCoef[,2],
         x1=obsCoef[,1] - 1.96 * obsCoef[,2])

obszero <- obsCoef[,1] - 1.96 * obsCoef[,2] < 0 &
           obsCoef[,1] + 1.96 * obsCoef[,2] > 0

points(x=obsCoef[,1], 
       y=1:nrow(obsCoef) - offset,
       pch=21, bg=ifelse(!obszero, "black", "white"))

axis(side=2, at=1:nrow(obsCoef), hadj=1, padj=0.5,
     labels=plotLabels, las=1)

text(x=relative.axis.point(0.015, "x"),
     y=relative.axis.point(0.97, "y"),
labels=expression(bold("(A)")*" Turnover"), adj=0)
box()
close.screen(1)

screen(2)
par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0))

plot(x=NULL, y=NULL, xlim=c(-0.03,0.03), ylim=rev(c(0.25, length(plot.order))), yaxt="n", xlab="",
     ylab="", xaxt="n")

# background dividers

rect(xleft=par("usr")[1], xright=par("usr")[2],
     ybottom=bDivs[-length(bDivs)]+0.5, ytop=bDivs[-1]+0.5,
     col=c("grey90","white"), border=NA)

segments(x0=0, x1=0, y0=par("usr")[3], y1=relative.axis.point(0.95,"y"),
         lty="31")
axis(side=1, mgp=c(3,0.1,0))

dRrands <- randRcoef[,1,]
obsCoef <- dRcoef

# rearrange to right order
obsCoef <- obsCoef[match(plot.order, rownames(obsCoef)),]
dRrands <- dRrands[match(plot.order, rownames(dRrands)),]

obszero <- obsCoef[,1] - 1.96 * obsCoef[,2] < 0 &
            obsCoef[,1] + 1.96 * obsCoef[,2] > 0

greyMean <- apply(dRrands, 1, mean)
colLims <- 0.015
greyMean <- (greyMean - -colLims) / (colLims - -colLims)
summary(greyMean)

sapply(1:nrow(dRrands), function(n){
  rect(xleft=quantile(dRrands[n,], prob=0.025),
       xright=quantile(dRrands[n,], prob=0.975),
       ybottom=n - 0.5 * boxWidth - offset,
       ytop = n  + 0.5 * boxWidth - offset,
       col=rgb(persRamp(greyMean[n])/255))
})

segments(y0=1:nrow(obsCoef) + offset, y1=1:nrow(obsCoef) + offset,
         x0=obsCoef[,1] + 1.96 * obsCoef[,2],
         x1=obsCoef[,1] - 1.96 * obsCoef[,2])

points(x=obsCoef[,1], 
       y=1:nrow(obsCoef) + offset,
       pch=21, bg=ifelse(!obszero, "black", "white"))

axis(side=2, at=1:nrow(obsCoef), 
     labels=NA, 
     las=1)

mtext(side=1, line=1, text="Standardized effect size",
      at=par("usr")[1])
box()

text(x=relative.axis.point(0.015, "x"),
     y=relative.axis.point(0.97, "y"),
      labels=expression(bold("(B)")*" Persistence"), adj=0)

close.screen(2)

# segments(x0=c(0.4,0.4), x1=c(0.6,0.6), y0=c(0.5,0.6), y1=c(0.5,0.6))
# points(y=c(0.5,0.6), x=c(0.5,0.5), pch=21, bg=c("white", "black"))
# text(x=0.6, y=c(0.5,0.6), labels=c("Non sig.", "Significant"), pos=4, offset=0.25)
# rect(xleft=0.4, xright=0.6, ybottom=0.365, ytop=0.435)
# par(lheight=0.85)
# text(x=0.6, y=0.4-0.01, labels="95% quantile of\nnull expectations", pos=4, offset=0.25)

close.screen(all.screens=TRUE)
dev.off()

#           interaction figure ####

pdf(date.wrap("./plots/post-novel interactions 3way", ".pdf"), height=3.35, width=8, bg="white")
par(mar=c(0,0,0,0), ps=8, mgp=c(3,0.5,0), tcl=-0.25, las=1)

ratScale <- c(-0.12, 0.18)

split.screen(rbind(c(0.1,0.35,0.35,0.95), # persistence N:dN
                   c(0.4,0.65,0.35,0.95), # persistence N:dpN
                   c(0.7,0.95,0.35,0.95), # persistence dN:dpN
                   
                   c(0.6,0.9,0.12,0.16), # persistence legend
                   
                   c(0.1,0.5,0.075,0.2))) # time series legend 
 
# novel observed distributions
screen(1)
par(mar=c(0,0,0,0), ps=8, mgp=c(3,0.5,0), tcl=-0.25, las=1)

richMat <- interaction.matrix.fun(dRdivDirDir, x.var="novelH1LS", y.var="deltaH1S", grid.size=100, quantile=0.025)

image(richMat, col=colorRampPalette(c("purple","white","orange"))(100), useRaster=TRUE,
      axes=FALSE, zlim=ratScale)
points(post200$novelH1LS ~ post200$deltaH1S, cex=0.2, pch=16, col="grey50")
contour(richMat, add=TRUE)
box()

richScale <- scaleDf[1:2,]

# delta axis
axis(side=2, at=(pretty(seq(-15,15,1),15)-richScale[2,1]) / richScale[2,2], 
     labels=pretty(seq(-15,15,1),15))
mtext(side=2, line=1.1, text = expression(Delta*"Diversity (pre-novel"%->%"novel)"),
      las=0)

# raw axis
rich <- seq(2,15,2)
axis(side=1, at=(rich-richScale[1,1]) / richScale[1,2], 
     labels=rich, mgp=c(3,0.1,0))
axis(side=1, at=(1:30-richScale[1,1]) / richScale[1,2], 
     labels=NA, tcl=-0.125)
mtext(side=1, line=0.75, text="Novel state diversity")

mtext(at=relative.axis.point(0, "x"), side=3, text=expression(bold("(A)")), adj=0)
close.screen(1)

screen(2)
par(mar=c(0,0,0,0), ps=8, mgp=c(3,0.5,0), tcl=-0.25, las=1)

richMat1 <- interaction.matrix.fun(dRdivDirDir, x.var="novelH1LS", y.var="deltaPostS", grid.size=100, quantile=0.025)

image(richMat1, col=colorRampPalette(c("purple","white","orange"))(100), useRaster=TRUE,
      axes=FALSE, zlim=ratScale)
points(post200$novelH1LS ~ post200$deltaPostS, cex=0.2, pch=16, col="grey50")
contour(richMat1, add=TRUE)
box()

richScale <- scaleDf[c(1,3),]

# delta axis
axis(side=2, at=(pretty(seq(-15,15,1),15)-richScale[2,1]) / richScale[2,2], 
     labels=pretty(seq(-15,15,1),15))
mtext(side=2, line=1.1, text = expression(Delta*"Diversity (novel"%->%"post-novel)"),las=0)

# raw axis
rich <- seq(2,15,2)
axis(side=1, at=(rich-richScale[1,1]) / richScale[1,2], 
     labels=rich, mgp=c(3,0.1,0))
axis(side=1, at=(1:30-richScale[1,1]) / richScale[1,2], 
     labels=NA, tcl=-0.125)
mtext(side=1, line=0.75, text="Novel state diversity")

mtext(at=relative.axis.point(0, "x"), side=3, text=expression(bold("(B)")), adj=0)

close.screen(2)

screen(3)
richMat2 <- interaction.matrix.fun(dRdivDirDir, x.var="deltaH1S", y.var="deltaPostS", grid.size=100, quantile=0.025)

image(richMat2, col=colorRampPalette(c("purple","white","orange"))(100), useRaster=TRUE,
      axes=FALSE, zlim=ratScale)
points(post200$deltaH1S ~ post200$deltaPostS, cex=0.2, pch=16, col="grey50")
contour(richMat2, add=TRUE)
box()

richScale <- scaleDf[2:3,]

# delta axis
axis(side=1, at=(pretty(seq(-15,15,1),15)-richScale[1,1]) / richScale[1,2], 
     labels=pretty(seq(-15,15,1),15), mgp=c(3,0.1,0))
mtext(side=1, line=0.75, text= expression(Delta*"Diversity (pre-novel"%->%"novel)"))
      
# raw axis
axis(side=2, at=(pretty(seq(-15,15,1),15)-richScale[2,1]) / richScale[2,2], 
     labels=pretty(seq(-15,15,1),15))

mtext(side=2, line=1.1, text= expression(Delta*"Diversity (novel"%->%"post-novel)"),las=0)

mtext(at=relative.axis.point(0, "x"), side=3, text=expression(bold("(C)")), adj=0)

close.screen(3)

screen(4)

image(y=c(0,1),
      x=seq(min(ratScale), max(ratScale), len=200),
      z=matrix(1:200, ncol=1),
      col=colorRampPalette(c("purple","white","orange"))(200),
      useRaster=TRUE, axes=FALSE, xlab="", ylab="")
axis(side=1, mgp=c(3,0,0))
mtext(side=1, text="Persistence", las=0, line=0.75)
box()
close.screen(4)

screen(5)
par(mar=c(0,0,0,0), ps=8, mgp=c(3,0.5,0), tcl=-0.25, las=1)
plot(NULL, xlim=c(0.25,0.75), ylim=c(0.2,0.8), axes=FALSE, xlab="", ylab="")

par(xpd=NA)
pointLocs = relative.axis.point(seq(0.15,0.85,len=3), "x")

Arrows(x0=c(relative.axis.point(0, "x"), pointLocs) + 0.03,
       x1=c(pointLocs, relative.axis.point(1, "x")) - 0.03,
       y0=0.5, y1=0.5, arr.type="triangle", arr.length = 0.15, arr.width=0.15)

points(x=pointLocs, y=rep(0.5,3),
       pch=c(24,21,21), bg=c("purple", "orange", "white"), cex=3)

segments(x0=pointLocs[2], x1=pointLocs[2], y0=0.6, y1=0.65)
text(x=pointLocs[2], y=0.65, pos=3, labels="Novel state diversity",offset=0.25)

segments(x0=pointLocs[c(1,3)], x1=pointLocs[c(1,3)], y0=0.4, y1=0.35)
segments(x0=pointLocs[2] + c(-0.01, 0.01), x1=pointLocs[2] + c(-0.01, 0.01), y0=0.4, y1=0.35)
segments(x0=pointLocs[c(1,3)], x1=pointLocs[2] + c(-0.01, 0.01), y0=0.35, y1=0.35)

text(x=mean(c(pointLocs[1], pointLocs[2] - 0.01)), 
     y=0.35, pos=1, adj=0.5, labels="Pre-novel to novel\ndiversity change")

text(x=mean(c(pointLocs[3], pointLocs[2] - 0.01)), 
     y=0.35, pos=1, adj=0.5, labels="Novel to post-novel\ndiversity change")
par(xpd=FALSE)
close.screen(5)

close.screen(all.screens=TRUE)

dev.off()

#           combo plot ####

library(viridisLite)

# dM rand coefs

pdf("./plots/turnPersEffects.pdf", height=5, width=10, useDingbats = FALSE)

par(mfrow=c(1,2), mar=c(3,3.5,1,0.5), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)

dMrands <- randMcoef[,1,][-(1:8),]
obsCoef <- dMcoef[-(1:8),]

coefCats = c(3,2,2,1,3,3,3,1,1,2)

plotLim = 0.135
plot(NULL, xlim=c(-0.05,0.15), ylim=c(-0.09,0.075),
     xlab="", ylab="", xaxt="n")
axis(side=1, mgp=c(3,0.1,0))

abline(h=0, col="grey")
abline(v=0,col="grey")
abline(a=0,b=1, lty="31")

text(x=relative.axis.point(0.165, "x"), y=0, col="grey", adj=0, pos=3, offset=0.25,
     labels="No post-novel effect")
text(y=relative.axis.point(0.975, "y"), x=0, col="grey", adj=0, pos=2, offset=0.25,
     labels="No null effect", srt=90)
text(x=relative.axis.point(0.125, "x"), y=relative.axis.point(0.125, "y"), col="black", adj=0, pos=3, offset=0.25,
     labels="Post-novel = Null", srt=45)

segments(x0 = rowMeans(dMrands), x1=rowMeans(dMrands),
         y0 = obsCoef[,1] + 1.96 * obsCoef[,2],
         y1 = obsCoef[,1] - 1.96 * obsCoef[,2], lwd=2,
         col=c("grey80","black","red")[coefCats])

segments(x0 = apply(dMrands, 1, quantile, probs=0.025), x1 = apply(dMrands, 1, quantile, probs=0.975),
         y0 = obsCoef[,1], y1 = obsCoef[,1], lwd=2,
         col=c("grey80","black","red")[coefCats])

points(obsCoef[,1] ~ rowMeans(dMrands),
       pch=16, col=c("grey80","black","red")[coefCats])

text(obsCoef[,1] ~ rowMeans(dMrands), labels=rownames(dMrands), pos=4,
     col=c("grey80","black","red")[coefCats])
mtext(side=1, line=1.25, text="Effect on turnover (non-novel null expectations)")
mtext(side=2, line=2.5, text="Effect on turnover (novel)", las=0)

text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.95, "y"),
     adj=0, labels=expression(bold("(A) ")*"Turnover"))



# dM rand coefs
dMrands <- randRcoef[,1,][-(1:8),]
obsCoef <- dRcoef[-(1:8),]

coefCats = c(3,2,3,1,2,2,3,2,3,3)

plotLim = 0.05
plot(NULL, xlim=c(-0.008,0.018), ylim=c(-0.02,0.025),
     xlab="", ylab="", xaxt="n")
axis(side=1, mgp=c(3,0.1,0))

abline(h=0, col="grey")
abline(v=0,col="grey")
abline(a=0,b=1, lty="31")

text(x=relative.axis.point(0.165, "x"), y=0, col="grey", adj=0, pos=3, offset=0.25,
     labels="No post-novel effect")
text(y=relative.axis.point(0.975, "y"), x=0, col="grey", adj=0, pos=2, offset=0.25,
     labels="No null effect", srt=90)
text(x=relative.axis.point(0.865, "x"), y=relative.axis.point(0.565, "y"), col="black", adj=1, pos=3, offset=0.25,
     labels="Post-novel = Null", srt=18)

segments(x0 = rowMeans(dMrands), x1=rowMeans(dMrands),
         y0 = obsCoef[,1] + 1.96 * obsCoef[,2],
         y1 = obsCoef[,1] - 1.96 * obsCoef[,2], lwd=2,
         col=c("grey80","black","red")[coefCats])

segments(x0 = apply(dMrands, 1, quantile, probs=0.025), x1 = apply(dMrands, 1, quantile, probs=0.975),
         y0 = obsCoef[,1], y1 = obsCoef[,1], lwd=2,
         col=c("grey80","black","red")[coefCats])

points(obsCoef[,1] ~ rowMeans(dMrands),
       pch=16, col=c("grey80","black","red")[coefCats])

text(obsCoef[,1] ~ rowMeans(dMrands), labels=rownames(dMrands), pos=4,
     col=c("grey80","black","red")[coefCats])
mtext(side=1, line=1.25, text="Effect on persistence (non-novel null expectations)")
mtext(side=2, line=2.5, text="Effect on persistence (novel)", las=0)

legend(x = relative.axis.point(0.625, "x"),
       y = relative.axis.point(0.15,"y"),
       pch = 22, pt.bg=c("grey80", "black", "red"),
       legend = c("No effect", "Post-novel = Null",
                  expression("Post-novel "!=" Null")),
       bty="n", pt.cex=1.5, y.intersp=0.7, x.intersp=0.85)

text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.95, "y"),
     adj=0, labels=expression(bold("(B) ")*"Persistence"))

dev.off()


#           combo plot (covariates) ####

library(viridisLite)

# dM rand coefs

pdf("./plots/turnPersEffectsCOVARIATES.pdf", height=5, width=10, useDingbats = FALSE)

par(mfrow=c(1,2), mar=c(3,3.5,1,0.5), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)

dMrands <- randMcoef[,1,][(2:8),]
obsCoef <- dMcoef[(2:8),]

coefCats = c(3,1,3,3,2,3,1)

plotLim = 0.135
plot(NULL, xlim=c(-0.1,0.1), ylim=c(-0.08,0.1),
     xlab="", ylab="", xaxt="n")
axis(side=1, mgp=c(3,0.1,0))

abline(h=0, col="grey")
abline(v=0,col="grey")
abline(a=0,b=1, lty="31")

text(x=relative.axis.point(0.165, "x"), y=0, col="grey", adj=0, pos=3, offset=0.25,
     labels="No post-novel effect")
text(y=relative.axis.point(0.975, "y"), x=0, col="grey", adj=0, pos=2, offset=0.25,
     labels="No null effect", srt=90)
text(x=relative.axis.point(0.125, "x"), y=relative.axis.point(0.125, "y"), col="black", adj=0, pos=3, offset=0.25,
     labels="Post-novel = Null", srt=45)

segments(x0 = rowMeans(dMrands), x1=rowMeans(dMrands),
         y0 = obsCoef[,1] + 1.96 * obsCoef[,2],
         y1 = obsCoef[,1] - 1.96 * obsCoef[,2], lwd=2,
         col=c("grey80","black","red")[coefCats])

segments(x0 = apply(dMrands, 1, quantile, probs=0.025), x1 = apply(dMrands, 1, quantile, probs=0.975),
         y0 = obsCoef[,1], y1 = obsCoef[,1], lwd=2,
         col=c("grey80","black","red")[coefCats])

points(obsCoef[,1] ~ rowMeans(dMrands),
       pch=16, col=c("grey80","black","red")[coefCats])

text(obsCoef[,1] ~ rowMeans(dMrands), labels=rownames(dMrands), pos=4,
     col=c("grey80","black","red")[coefCats])
mtext(side=1, line=1.25, text="Effect on turnover (non-novel null expectations)")
mtext(side=2, line=2.5, text="Effect on turnover (novel)", las=0)

text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.95, "y"),
     adj=0, labels=expression(bold("(A) ")*"Turnover"))

# dM rand coefs
dMrands <- randRcoef[,1,][(2:8),]
obsCoef <- dRcoef[(2:8),]

coefCats = c(3,1,1,3,1,2,1)

plotLim = 0.05
plot(NULL, xlim=c(-0.005,0.008), ylim=c(-0.013,0.03),
     xlab="", ylab="", xaxt="n")
axis(side=1, mgp=c(3,0.1,0))

abline(h=0, col="grey")
abline(v=0,col="grey")
abline(a=0,b=1, lty="31")

text(x=relative.axis.point(0.165, "x"), y=0, col="grey", adj=0, pos=3, offset=0.25,
     labels="No post-novel effect")
text(y=relative.axis.point(0.975, "y"), x=0, col="grey", adj=0, pos=2, offset=0.25,
     labels="No null effect", srt=90)
text(x=relative.axis.point(0.865, "x"), y=relative.axis.point(0.565, "y"), col="black", adj=1, pos=3, offset=0.25,
     labels="Post-novel = Null", srt=18)

segments(x0 = rowMeans(dMrands), x1=rowMeans(dMrands),
         y0 = obsCoef[,1] + 1.96 * obsCoef[,2],
         y1 = obsCoef[,1] - 1.96 * obsCoef[,2], lwd=2,
         col=c("grey80","black","red")[coefCats])

segments(x0 = apply(dMrands, 1, quantile, probs=0.025), x1 = apply(dMrands, 1, quantile, probs=0.975),
         y0 = obsCoef[,1], y1 = obsCoef[,1], lwd=2,
         col=c("grey80","black","red")[coefCats])

points(obsCoef[,1] ~ rowMeans(dMrands),
       pch=16, col=c("grey80","black","red")[coefCats])

text(obsCoef[,1] ~ rowMeans(dMrands), labels=rownames(dMrands), pos=4,
     col=c("grey80","black","red")[coefCats])
mtext(side=1, line=1.25, text="Effect on persistence (non-novel null expectations)")
mtext(side=2, line=2.5, text="Effect on persistence (novel)", las=0)

legend(x = relative.axis.point(0.625, "x"),
       y = relative.axis.point(0.15,"y"),
       pch = 22, pt.bg=c("grey80", "black", "red"),
       legend = c("No effect", "Post-novel = Null",
                  expression("Post-novel "!=" Null")),
       bty="n", pt.cex=1.5, y.intersp=0.7, x.intersp=0.85)

text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.95, "y"),
     adj=0, labels=expression(bold("(B) ")*"Persistence"))

dev.off()

# -------------- ####
# POST-HOC TESTS ####
# -------------- ####
#       Exploring random intercepts ####
dMri = ranef(dMdivAbsAbs)[[1]]
dMri = cbind(dMri, tapply(predict(dMdivAbsAbs), post200$novelID, mean))
dMri = cbind(dMri, tapply(resid(dMdivAbsAbs), post200$novelID, mean))
dMri = cbind(dMri, tapply(post200$dMag, post200$novelID, mean))
colnames(dMri) = c("dMri", "dMpred", "dMres", "dM")
dMri$novelID = rownames(dMri)

dRri = ranef(dRdivDirDir)[[1]]
dRri = cbind(dRri, tapply(predict(dRdivDirDir), post200$novelID, mean))
dRri = cbind(dRri, tapply(resid(dRdivDirDir), post200$novelID, mean))
dRri = cbind(dRri, tapply(post200$dRatio, post200$novelID, mean))
colnames(dRri) = c("dRri", "dRpred", "dRres", "dR")
dRri$novelID = rownames(dRri)

cor(dMri[,(1:4)])
cor(dRri[,(1:4)])

plot((dMri$dMri + dMri$dMres + dMri$dMpred) ~ dMri$dM)
plot((dRri$dRri + dRri$dRres + dRri$dRpred) ~ dRri$dR)

post200 = merge(post200, dMri,
                by.x="novelID", by.y="novelID",
                all.x=TRUE, all.y=FALSE, sort=FALSE)

post200 = merge(post200, dRri,
                by.x="novelID", by.y="novelID",
                all.x=TRUE, all.y=FALSE, sort=FALSE)

pdf("./plots/randomInterceptCor.pdf", height=4.5, width=5, useDingbats = FALSE)
par(mar=c(2.5,3,0.5,0.5), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(dRri$dRri ~ dMri$dMri, type="n", xaxt="n")

rilm <- lm(dRri$dRri ~ dMri$dMri)
riPred = predict(rilm, se.fit=TRUE)

axis(side=1, mgp=c(3,0.1,0))
mtext(side=1, line=1, text="Turnover random intercept")
mtext(side=2, line=2, text="Persistence random intercept", las=0)

points(dRri$dRri ~ dMri$dMri, pch=16, cex=0.5, col="grey")

lines(x=sort(dMri$dMri), y=riPred$fit[order(dMri$dMri)], lwd=2)
polygon(x=c(sort(dMri$dMri), rev(sort(dMri$dMri))), 
         y=c(riPred$fit[order(dMri$dMri)] + 1.96 * riPred$se.fit[order(dMri$dMri)],
             rev(riPred$fit[order(dMri$dMri)] - 1.96 * riPred$se.fit[order(dMri$dMri)])),
         border=NA, col=rgb(0.5,0.5,0.5,0.25))

abline(h=0, col="grey", lty="31")
abline(v=0, col="grey", lty="31")
dev.off()

cor(cbind(dRri$dRri,dMri$dMri))

dirTable = table(dRri$dRri > 0, dMri$dMri > 0)
dirTable / sum(dirTable)
       
# with random models

randomRIs = do.call("rbind", lapply(rNgamList, function(x){
  
Ris <- cbind(ranef(x$Mfull)[[1]],
             ranef(x$Rfull)[[1]])

mR2diff = performance(x$Mfull)
mR2diff = mR2diff$R2_conditional - mR2diff$R2_marginal

rR2diff = performance(x$Rfull)
rR2diff = rR2diff$R2_conditional - rR2diff$R2_marginal

return(c(mR2diff = mR2diff,
         rR2diff = rR2diff,
         riCor = cor(Ris)[2,1],
         riSlope = summary(lm(Ris[,2] ~ Ris[,1]))$coefficients[2,1],
         riR2 = performance(lm(Ris[,2] ~ Ris[,1]))$R2))
}))

randomMeans = apply(randomRIs,2, mean)
randomQuants = apply(randomRIs,2, quantile, probs=c(0.025,0.975))

diff(unlist(performance(dMdivAbsAbs)[5:4])) / randomMeans[1]
diff(unlist(performance(dRdivDirDir)[5:4])) / randomMeans[2]

#       Global vs regional temp change ####

tempComb = merge(traceLong, temp.df,
                 by.x="bin", by.y="age",
                 all.x=TRUE, all.y=FALSE, sort=FALSE)

tempSamp = tempComb[sample(1:nrow(tempComb), 1e5),]
plot(tempSamp$localTemp ~ tempSamp$fit)

# best-fitting change

library(viridisLite)

with(post200[!is.na(post200$lat),],
plot(localLag25 ~ globalLag25, pch=21,
     bg=rgb(colorRamp(viridis(5))(unitScale(lat))/255)))
abline(h=0)
abline(v=0)
with(post200[!is.na(post200$lat),],
     text(localLag25 ~ globalLag25, labels=round(lat), pos=4, cex=0.6))

cor(post200[,c('localLag5','localLag10','localLag15','localLag20', 'localLag25', 
               'globalLag5','globalLag10','globalLag15','globalLag20','globalLag25')])

#       Extract site differences ####


siteDescr = site.df$description
hasPhysio = grepl("Physiography:", siteDescr)
hasVeg = grepl("Surrounding vegetation:|Vegetation formation:", siteDescr)

hasStuff = siteDescr[hasPhysio & hasVeg]
hasStuff = gsub("Vegetation formation:", "Surrounding vegetation:", hasStuff)

sitePhys = substr(hasStuff, regexpr("Physiography:", hasStuff) + nchar("Physiography:") + 1, regexpr("Surrounding vegetation:", hasStuff)-3)
sort(table(tolower(sitePhys)), decreasing=TRUE)


#       gam of turnover and persistence axes ####

# fixed & random residuals
post200$dMagRes = residuals(dMdivAbsAbs)
post200$dRatioRes = residuals(dRdivDirDir)

# fixed only residuals
post200$dMagResFixed = log(post200$dMagRaw) - predict(dMdivAbsAbs, re.form=NA)
post200$dRatioResFixed = post200$dRatio - predict(dRdivDirDir, re.form=NA)

# model predictions (fixed only)
post200$dMagPredFixed = predict(dMdivAbsAbs, re.form=NA)
post200$dRatioPredFixed = predict(dRdivDirDir, re.form=NA)

# model predictions (fixed + random)
post200$dMagPred = predict(dMdivAbsAbs)
post200$dRatioPred = predict(dRdivDirDir)

# random intercepts
riDf = data.frame(magRI = ranef(dMdivAbsAbs)[[1]][,1],
                  ratioRI = ranef(dRdivDirDir)[[1]][,1])

pdf("./plots/modelPreds.pdf", height=5, width=6, useDingbats=FALSE)
par(mfcol=c(2,2), mar=c(2,3,0.5,0.5), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)

sapply(1:4, function(n){
  
  if(n != 2){
  magVar = post200[,c("dMagResFixed", "dMagPred","dMagPredFixed", "dMagPred")][,n]
  ratVar = post200[,c("dRatioResFixed", "dRatioPred","dRatioPredFixed", "dRatioPred")][,n]
  } else {
  magVar = riDf[,1]
  ratVar = riDf[,2]
  }
  tempGam = gam(ratVar ~ s(magVar))
  
  tempPred = data.frame(magVar = seq(min(magVar), max(magVar), len=200))
  tempPred = cbind(tempPred,
                   as.data.frame(predict(tempGam, newdata=tempPred, se.fit=TRUE)))
  
  if(n %in% c(3:4)){
    magVar = exp(magVar)
    tempPred$magVar = exp(tempPred$magVar)
    }

  
  plot(ratVar ~ magVar, pch=16, cex=0.5, col="grey85",
       ylim=quantile(ratVar, probs=c(0.01, 0.99)), xaxt="n", xlab="", ylab="")
  
  axis(side=1, mgp=c(3,0.1,0))
  
  if(n %in% c(1,2)){
    mtext(side=1, line=1, text="Residual turnover", cex=0.8)
    mtext(side=2, line=2, text="Residual persistence", las=0, cex=0.8)
  } else {
    mtext(side=1, line=1, text="Turnover", cex=0.8)
    mtext(side=2, line=2, text="Persistence", las=0, cex=0.8)
  }
  
  abline(h=0, lty="31")
  
  rampRange = max(abs(range(tempPred$fit)))
  
  tempPred$predRamp = unitScale(tempPred$fit, custMin= -rampRange, custMax= rampRange)
  
  sapply(2:nrow(tempPred), function(n){
    
    rampCol = colorRamp(c("purple","grey70","orange"))(tempPred$predRamp[n])/255
    
    lines(tempPred$fit[(n-1):n] ~ tempPred$magVar[(n-1):n], lwd=2,
          col=rgb(rampCol[1], rampCol[2], rampCol[3]))
    
    polygon(x=c(tempPred$magVar[(n-1):n], rev(tempPred$magVar[(n-1):n])),
            y=c(tempPred$fit[(n-1):n] + 1.96 * tempPred$se.fit[(n-1):n],
                rev(tempPred$fit[(n-1):n] - 1.96 * tempPred$se.fit[(n-1):n])),
            col=rgb(rampCol[1], rampCol[2], rampCol[3], 0.5),
            border=NA)
    
    
  })
  
  text(x=relative.axis.point(0.025, "x"),
       y=relative.axis.point(0.95, "y"),
       labels=paste0("(", LETTERS[c(1,3,2,4)][n], ")"), 
       font=2, adj=0)
  
  text(x=relative.axis.point(0.105, "x"),
       y=relative.axis.point(0.95, "y"),
       labels=c("Fixed effect residuals",
                "Novel state random intercepts",
                "Fixed effect predictions",
                "Full model predictions")[n],
       font=1, adj=0)
  
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.05, "y"),
       labels=expression("R"^2*" = "), adj=0)
  
  text(x=relative.axis.point(0.13, "x"),
       y=relative.axis.point(0.04, "y"),
       labels=sprintf("%.3f", performance(tempGam)$R2), adj=0)
  
    
})
dev.off()

# random intercepts
riDf = data.frame(magRI = ranef(dMdivAbsAbs)[[1]][,1],
                  ratioRI = ranef(dRdivDirDir)[[1]][,1])

# models
MRresGam = gam(dRatioResFixed ~ s(dMagResFixed) + time.since.novelS, data=post200)
riGam = gam(ratioRI ~ s(magRI), data=riDf)

par(mfrow=c(1,3))
# raw plot
plot(post200$dRatio ~ post200$dMag, pch=16, col="grey80")
abline(h=0, lty="31")
ratPred = seq(min(post200$dMag), max(post200$dMag), len=200)
pred.df = cbind(ratPred,
                as.data.frame(predict(MRgam, newdata=data.frame(dMag=ratPred, time.since.novelS=0), se.fit=TRUE)))
polygon(x=c(pred.df$ratPred, rev(pred.df$ratPred)),
        y=c(pred.df$fit + 1.96 * pred.df$se.fit,
            rev(pred.df$fit - 1.96 * pred.df$se.fit)),
        col=rgb(1,0,0,0.5), border=NA)
lines(pred.df$fit ~ pred.df$ratPred, col="red", lwd=2)                

# fixed res plot 
plot(post200$dRatioResFixed ~ post200$dMagResFixed, pch=16, col="grey80")
abline(h=0, lty="31")
abline(v=0, lty="31")
ratPred = seq(min(post200$dMagResFixed), max(post200$dMagResFixed), len=200)
pred.df = cbind(ratPred,
                as.data.frame(predict(MRresGam, newdata=data.frame(dMagResFixed=ratPred, time.since.novelS=0), se.fit=TRUE)))
polygon(x=c(pred.df$ratPred, rev(pred.df$ratPred)),
        y=c(pred.df$fit + 1.96 * pred.df$se.fit,
            rev(pred.df$fit - 1.96 * pred.df$se.fit)),
        col=rgb(1,0,0,0.5), border=NA)
lines(pred.df$fit ~ pred.df$ratPred, col="red", lwd=2)

# RI plot
plot(riDf$ratioRI ~ riDf$magRI, col="grey80", pch=16)
abline(h=0, lty="31")
abline(v=0, lty="31")
ratPred = seq(min(riDf$magRI), max(riDf$magRI), len=200)
pred.df = cbind(ratPred,
                as.data.frame(predict(riGam, newdata=data.frame(magRI=ratPred), se.fit=TRUE)))
polygon(x=c(pred.df$ratPred, rev(pred.df$ratPred)),
        y=c(pred.df$fit + 1.96 * pred.df$se.fit,
            rev(pred.df$fit - 1.96 * pred.df$se.fit)),
        col=rgb(1,0,0,0.5), border=NA)
lines(pred.df$fit ~ pred.df$ratPred, col="red", lwd=2)


#       holocene models (restrict novel state occurrence timeframe) ####

holoData <- post200[post200$novel.time <= 9700,]

holoMagM <- update(dMdivAbsAbs, .~., data=holoData)
holoRatioM <- update(dRdivDirDir, .~., data=holoData)

holoMagCI = cbind(summary(holoMagM)$coefficients[,1],
                  summary(holoMagM)$coefficients[,1] - 1.96 * summary(holoMagM)$coefficients[,2],
                  summary(holoMagM)$coefficients[,1]+ 1.96 * summary(holoMagM)$coefficients[,2])
allMagCI = cbind(summary(dMdivAbsAbs)$coefficients[,1],
                  summary(dMdivAbsAbs)$coefficients[,1] - 1.96 * summary(dMdivAbsAbs)$coefficients[,2],
                  summary(dMdivAbsAbs)$coefficients[,1]+ 1.96 * summary(dMdivAbsAbs)$coefficients[,2])

holoRatioCI = cbind(summary(holoRatioM)$coefficients[,1],
                  summary(holoRatioM)$coefficients[,1] - 1.96 * summary(holoRatioM)$coefficients[,2],
                  summary(holoRatioM)$coefficients[,1]+ 1.96 * summary(holoRatioM)$coefficients[,2])
allRatioCI = cbind(summary(dRdivDirDir)$coefficients[,1],
                 summary(dRdivDirDir)$coefficients[,1] - 1.96 * summary(dRdivDirDir)$coefficients[,2],
                 summary(dRdivDirDir)$coefficients[,1]+ 1.96 * summary(dRdivDirDir)$coefficients[,2])

pdf("./plots/holoceneCoefComp.pdf", height=4, width=8, useDingbats = FALSE)

par(mfrow=c(1,2), mar=c(0,2.5,0,0), oma=c(3,1.25,1,1), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)

plot(holoMagCI[,1] ~ allMagCI[,1],
     xlim=c(-0.1,0.1), ylim=c(-0.11,0.1), type="n", ylab="", xlab="", axes=FALSE)
axis(side=1, mgp=c(3,0.2,0))
axis(side=2)
mtext(side=2, line=2.5, text="Fixed effect coefficients (3,000 to 9,700 ybp)", las=0)

abline(a=0,b=1, lty="dashed")

segments(x0=allMagCI[,1], x1=allMagCI[,1],
         y0=holoMagCI[,3], y1=holoMagCI[,2],
         col=c(rep("grey80", 10), rep("black", 2), rep("grey80", 8)))
segments(x0=allMagCI[,2], x1=allMagCI[,3],
         y0=holoMagCI[,1], y1=holoMagCI[,1],
         col=c(rep("grey80", 10), rep("black", 2), rep("grey80", 8)))
points(holoMagCI[,1] ~ allMagCI[,1], pch=16,
       col=c(rep("grey80", 10), rep("black", 2), rep("grey80", 8)))

text((holoMagCI[11:12,1] + 0.005) ~ allMagCI[c(11:12),1],
     pos=4, labels=c("Global temperature change", "Regional temperature change"))

text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.95, "y"),
     adj=0, labels=expression(bold("(A)")*" Post-novel turnover"), font=2)

box()

plot(holoRatioCI[,1] ~ allRatioCI[,1],
     xlim=c(-0.02,0.04), ylim=c(-0.02,0.04), type="n", ylab="", xlab="", axes=FALSE)
axis(side=1, mgp=c(3,0.2,0))
axis(side=2)

abline(a=0,b=1, lty="dashed")

mtext(side=1, line=1, text="Fixed effect coefficients (3,000 to 19,000 ybp)", at=par("usr")[1])
segments(x0=allRatioCI[,1], x1=allRatioCI[,1],
         y0=holoRatioCI[,3], y1=holoRatioCI[,2],
         col=c(rep("grey80", 10), rep("black", 2), rep("grey80", 8)))
segments(x0=allRatioCI[,2], x1=allRatioCI[,3],
         y0=holoRatioCI[,1], y1=holoRatioCI[,1],
         col=c(rep("grey80", 10), rep("black", 2), rep("grey80", 8)))
points(holoRatioCI[,1] ~ allRatioCI[,1], pch=16,
       col=c(rep("grey80", 10), rep("black", 2), rep("grey80", 8)))

text((holoRatioCI[11:12,1] + 0.0015) ~ allRatioCI[c(11:12),1],
     pos=4, labels=c("Global temperature change", "Regional temperature change"))

text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.95, "y"),
     adj=0, labels=expression(bold("(B)")*" Post-novel persistence"), font=2)

box()
dev.off()

#       climate lag model comparison ####

climAIC$gLag <- as.numeric(gsub("globalLag|S","",climAIC$globalLag))
climAIC$lLag <- as.numeric(gsub("localLag|S","",climAIC$localLag))

library(viridisLite)
pdf("./plots/climateLagAIC.pdf", height=4.5, width=4.5, useDingbats = FALSE)

par(mar=c(3,4,0.5,0.5), ps=10, tcl=-0.25, 
    mgp=c(3,0.5,0), las=1)

plot(climAIC$dR ~ climAIC$dM, type="n", xlab="", ylab="", xaxt="n",
     xlim=range(climAIC$dM) + c(0,5))
axis(side=1, mgp=c(3,0.2,0))
mtext(side=1, line=1.5, text=expression("AIC (post-novel turnover)"))
mtext(side=2, line=3, text=expression("AIC (post-novel persistence)"), las=0)

sapply(split(climAIC, f=climAIC$gLag),
       function(x){
      
         lines(x$dR ~ x$dM)
         points(x$dR ~ x$dM,
                bg=viridis(6)[as.factor(x$lLag)], pch=21)
            
         text(x=x$dM[x$lLag==1],
              y=x$dR[x$lLag==1],
              pos=4,
              labels=paste0("Gl = ", format(x$gLag[x$lLag==1]*200, big.mark=",")))
         
         if(x$gLag[1]==10){
           
           text(x=mean(x$dM),
                y=max(x$dR),
                pos=3, offset=1, labels="Local lag")
           
           text(x=x$dM, y=x$dR,
                pos=3, labels=format(x$lLag*200, big.mark=","))
           
           
         }
       })

dev.off()

# ------------- ####
# SUMMARY STATS ####

# make predictions for time from novel at 200-1000 years
mainCoef <- rownames(obsCoef)[!grepl(":", rownames(obsCoef))]
mPred <- as.data.frame(matrix(0, ncol=length(mainCoef), nrow=5,
                              dimnames=list(seq(200,1000,200), mainCoef)))
mPred[,"time.since.novelS"] = (seq(200,1000,200) - scaleDf["time.since.novel",1]) / scaleDf["time.since.novel", 2]

mPred <- predict(dMdiv, newdata=mPred, re.form=NA)

table(sapply(split(novelTrajectories[[1]][[3]],
             f=novelTrajectories[[1]][[3]]$novelID), function(x){min(x$time.since.novel) <= 1000}))

templower <- postNovDiffs[[1]]$obs
templower[upper.tri(templower, diag=TRUE)] = NA
tempupper <- postNovDiffs[[1]]$obs
tempupper[lower.tri(tempupper, diag=TRUE)] = NA
tempdiag <- postNovDiffs[[1]]$obs
tempdiag[lower.tri(tempdiag, diag=FALSE)] = NA
tempdiag[upper.tri(tempdiag, diag=FALSE)] = NA

tempobs <- postNovDiffs[[1]]$obs
tempprop <- postNovDiffs[[1]]$propMean

# more total
sum(postNovDiffs[[1]]$obs[tempprop > 1 & !is.na(tempprop)])

sum(tempupper[tempprop > 1 & !is.na(tempprop)], na.rm=TRUE) # more upper
sum(templower[tempprop > 1 & !is.na(tempprop)], na.rm=TRUE) # more lower
sum(tempdiag[tempprop > 1 & !is.na(tempprop)], na.rm=TRUE) # more diag

# less total
sum(postNovDiffs[[1]]$obs[tempprop < 1 & !is.na(tempprop)])

sum(tempupper[tempprop < 1 & !is.na(tempprop)], na.rm=TRUE) # more upper
sum(templower[tempprop < 1 & !is.na(tempprop)], na.rm=TRUE) # more lower
sum(tempdiag[tempprop < 1 & !is.na(tempprop)], na.rm=TRUE) # more diag

temprand <- postNovDiffs[[1]]$randMean

sum(temprand[postNovDiffs[[1]]$obs == 0] > 0)

# novel persistence and temperature
plot(post200$dRatio ~ post200$tempLS)

# Neotoma plant obs
sub.plant <- droplevels(plant.record.df[plant.record.df$site.1 %in% post200$site,])
length(unique(sub.plant$))

novDf <- do.call("rbind", all.novel$novel)
sum(novDf$novel, na.rm=TRUE)
length(unique(novDf$site[novDf$novel]))
length(unique(post200$site))

plot(post200$gamma ~ post200$deltaH1)

x = unique(post200$site)
df = do.call("rbind", all.novel$novel)
df = df[complete.cases(df),]
sum(df$novel, na.rm=TRUE) ; length(unique(df$site)) ; nrow(df)

df = df[as.numeric(df$bins) >= 3000 & as.numeric(df$bins) < 17000,]
sum(df$novel, na.rm=TRUE) ; length(unique(df$site)) ; nrow(df)

sum(df$novel, na.rm=TRUE) - length(unique(post200$novelID))

a = table(post200$dRatio >0, post200$time.since.novel)
a[2,] / colSums(a)

# diversity divisions

plot((post200$deltaH1 - post200$deltaPost) ~ post200$H1,
     pch=16)
abline(h=0, lty="31", col="grey")

table(post200$deltaH1 < 0, post200$deltaPost < 0)

# do an RDA for which taxa load best with lower diversity?
matTax = sort(unique(unlist(sapply(all.novel$prop.ssmats, colnames))))

allNovMat = do.call('rbind', lapply(1:length(all.novel$prop.ssmats), function(n){
  # add on missing columns
  x = all.novel$prop.ssmats[[n]]
  templateMat = matrix(0, nrow=nrow(x), ncol=length(matTax), 
                       dimnames=list(paste0(all.novel$site[n], ":", rownames(x)), matTax))
  templateMat[,match(colnames(x), colnames(templateMat))] = x
  return(templateMat)
  }))

# now just the ones in our post200 dataset
post200Mat = allNovMat[match(paste0(post200$site, ":", post200$bin), rownames(allNovMat)),]

tempRDA = dbrda(post200Mat ~ H1, data=post200)
sppscores(tempRDA) <- sqrt(decostand(post200Mat, "total"))
plot(tempRDA, type="text")
summary(tempRDA)$species

# was diversity loss just loss of evenness?
post200$deltaH0 = post200$novelH0 - post200$preNovH0
post200$deltaPostH0 = post200$H0 - post200$novelH0

plot(post200$deltaH1 ~ post200$deltaH0)
table(post200$deltaH1 < 0, post200$deltaH0 < 0)

plot(post200$deltaPost ~ post200$deltaPostH0)
table(post200$deltaPost < 0, post200$deltaPostH0 < 0)

hist(post200$globalLag25)
table(post200$globalLag25 > 0)
length(unique(post200$bin[post200$globalLag25 <= 0]))

hist(post200$localLag15)
table(post200$localLag15 > 0)
length(unique(post200$bin[post200$localLag15 <= 0]))

summary(post200$novelH0[!duplicated(post200$novelID)])
sd(post200$novelH0)
summary(post200$deltaH1[!duplicated(post200$novelID)] < 0)

table(post200$deltaH1 > 0, post200$deltaPost > 0)
table(post200$deltaH1 > 0, post200$deltaPost > 0) / nrow(post200)

# ------------- ####
#           compare raw count and PPE-corrected novelty ####

ppeNovel = do.call("rbind", readRDS("./outputs/all neotoma novelty.rds")$novel)
ppeSub = ppeNovel[,c("site", "bins", "seq.dist", "raw.min.dist", "seq.p", "min.p", "cat")]
colnames(ppeSub)[-(1:2)] = paste0("ppe.", colnames(ppeSub)[-(1:2)])

rawNovel = do.call("rbind", readRDS("./outputs/all neotoma novelty.raw.rds")$novel)
rawNovel = merge(rawNovel, ppeSub, by.x=c("site", "bins"), by.y=c("site", "bins"),
                 all.x=TRUE, all.y=FALSE, sort=FALSE)

cor(rawNovel[,c("seq.dist", "ppe.seq.dist")], use="complete.obs")
cor(rawNovel[,c("raw.min.dist", "ppe.raw.min.dist")], use="complete.obs")
cor(rawNovel[,c("seq.p", "ppe.seq.p")], use="complete.obs")
cor(rawNovel[,c("min.p", "ppe.min.p")], use="complete.obs")

table(rawNovel$cat, rawNovel$ppe.cat)

#           temperature change correlation matrix ####

pdf("./plots/localGlobalTempComp.pdf", height=5, width=6.5, useDingbats = FALSE)
par(mar=c(3,3,0.5,0.5), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(post200$localLag25 ~ post200$globalLag25, pch=16,
     col=rgb(colorRamp(viridis(5))(unitScale(abs(post200$lat)))/255), xaxt="n")
abline(h=0, col="grey90")
axis(side=1, mgp=c(3,0.2,0))
mtext(side=1, line=1, text=expression("Global temperature change ("*Delta*degree*"C)"))
mtext(side=2, line=2, text=expression("Regional temperature change ("*Delta*degree*"C)"), las=0)
dev.off()




tempCols = colnames(post200)[grepl("Lag", colnames(post200)) & !grepl("S", colnames(post200))]
write.csv(cor(post200[,sort(tempCols)]), "./outputs/tempChangeCorMat.csv")

#           diversity correlation (pre, novel, post) ####

cor(post200[,c("novelH0", "novelH1", "novelH2", "preNovH0", "preNovH1", "preNovH2", "H0", "H1", "H2")])

#           diversity correlation (H0, H1, H2, rarefaction) ####

novTrajDf = novelTrajectories[[1]]$novel
cor(novTrajDf[,c("H0","H1","H2", "rareType","rareDiv")])

backTrajDf = novelTrajectories[[2]]
backTrajDf = backTrajDf[!duplicated(backTrajDf[,c("site", "bin")]),]

divCor = as.matrix(cor(backTrajDf[,c("H0","H1","H2", "rareType","rareDiv")]))
divCor = apply(round(divCor, 4), 1, function(x){sprintf("%.3f", x)})
write.csv(divCor, "./outputs/diversityCorMat.csv")

# what do you know? Doesn't matter what you do, diversity estimates are all
# very similar.
dev.off()
plot(backTrajDf$rareDiv ~ backTrajDf$H0, log="xy")
plot(backTrajDf$rareType ~ backTrajDf$H1, log="x")

