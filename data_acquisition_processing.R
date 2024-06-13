# ################################# ####
# Predicting Novel Communities      ####
# Author:    Timothy L Staples      ####
# Collaborators: John Pandolfi      #### 
#                Wolfgang Kiessling ####
# ################################# ####
# Script purpose ####
# Global attributes & working directories ####

rm(list=ls())
setwd("/Users/uqtstapl/Library/CloudStorage/Dropbox/Tim/Post-doc/Research projects/novel_comms_legacy/prodCode")

# Global functions ####

# source functions from 'functions' sub-folder
sapply(paste0("./functions/", list.files("./functions", pattern =".R")), source)

# Packages ####

package.loader(c("neotoma2", "WorldFlora"))

# A little custom function to add dates to output files
date.wrap <- function(string, ext){
  paste0(string, " ", Sys.Date(), ext)
}

# ####
# 1. DATA IMPORT <- This takes a while ####
#           Download from Neotoma ####

# Download all vascular plant data from Neotoma, then save output (nested list) as 
# and .rds file.

plant.datasets = get_datasets(datasettype = "pollen", all_data = TRUE)

plantRecords = lapply(sort(unique(c(seq(1,length(plant.datasets),10), length(plant.datasets))))[-1],
                      function(n){
                        print(paste0((n-10), " - ", n))
                        return(get_downloads(plant.datasets[(n-10):n]))
                      })

saveRDS(plantRecords, file="./rawdata/Neotoma vascular plant records.rds")

#           Convert database tables into complete data-frame ####

plantRecords <- readRDS("./rawdata/Neotoma vascular plant records.rds")

# Convert relational SQL-form data (as an S4 object) to a data-frame with all data
# replicated on each row.

plantSamples = do.call("rbind", lapply(plantRecords, function(x){samples(x)}))

# This is a big file, so I'm saving it to recall it to save having to redo this
# section of code each time
write.csv(plantSamples, "./rawdata/Neotoma vascular plant records.csv")

# ####
# 2. TAXONOMY CLEANING ####

#           Load Neotoma data ####

# read in imported data if it's not already in the global environment
if(! "plantSamples" %in% ls()){
plant.record.df <- read.csv("./rawdata/Neotoma vascular plant records.csv")
}

#           Clean taxonomy ####
#                               Create taxa reference table ####

# Get taxa master data-frame, with ID numbers so we can match our 
# final results back to our plant record dataframe
taxa.table = as.data.frame(table(plantSamples$variablename))
taxa.table = taxa.table[order(taxa.table$Freq, decreasing=TRUE),]
colnames(taxa.table) = c("taxa", "count")

write.csv(taxa.table, "./outputs/rawTaxaTable.csv")

#                               Create word-level taxa table ####

# By splitting each word in the taxon field, we can search for that specific
# word in lists of genera, families etc. This specific process lets us do this
# in a vectorized form which is much faster than scanning through words using
# substring sampling or using loops.

WFO.remember(WFO.file = "./rawdata/WFO_Backbone.zip", WFO.data = "WFO.data", WFO.pos = 1)

taxaProc = gsub("[[:punct:]]&&[^/]]|undiff|type|subg|\\.", "", taxa.table$taxa)
taxaProc = gsub("/|-", " ", taxaProc)
taxaProc = gsub("ë", "e", taxaProc)
taxa.words <- strsplit(taxaProc, " ")

# some rows with "/" don't work correctly, especially if there's been a 
# genus/family change (e.g., Juniperus is "accepted" in Cupresseaceae and Podocarpaceae)
# THIS WILL TAKE TIME TO RUN
taxMatch = do.call("rbind", lapply(taxa.words, function(w){
  
  print(paste0(w, collapse = " "))
  taxRow = data.frame(family=NA, genus=NA, species=NA)
  
  # test for a full match
  tempMatch = WFO.match(spec.data=paste0(w, collapse = " "), 
                        WFO.data = WFO.data, Fuzzy=0, Fuzzy.force=FALSE,
                        spec.name.nobrackets=TRUE, verbose=FALSE)
  
  # to start with, keep only accepted
  tempAccept = WFO.one(tempMatch)
  
  # if we got a perfect match
  if(tempAccept$Matched & !is.na(tempAccept$Matched)){
    taxRow[1,"family"] = tempAccept$family
    taxRow[1,"genus"] = tempAccept$genus
    taxRow[1,"species"] = tempAccept$specificEpithet
    return(taxRow)
  }
  
  # if we didn't get a proper match, and we have more than one word, we can
  # try for word-based searches. This will catch  etc.
  if(length(w) > 1){
    print("Matching each word")
    wordMatch = do.call("rbind", lapply(w, function(wS){
      tempMatch = WFO.match(spec.data=wS, 
                            WFO.data = WFO.data, Fuzzy=0, Fuzzy.force=FALSE,
                            spec.name.nobrackets=TRUE, verbose=FALSE)
    }))
    
    wordMatch = wordMatch[wordMatch$Matched,]
    
    # Only one word matches (E.g., "Lycopodium spike")
    if(sum(wordMatch$Matched) == 1){
      taxRow[1,"family"] = wordMatch$family[wordMatch$Matched]
      taxRow[1,"genus"] = wordMatch$genus[wordMatch$Matched]
      taxRow[1,"species"] = wordMatch$specificEpithet[wordMatch$Matched]
      return(taxRow)
    }
    
    # More than one word matches (E.g., "Ostrya/Carpinus")
    # We retain only the upper level details
    if(sum(wordMatch$Matched) > 1){
      taxRow[1,"family"] = ifelse(length(unique(wordMatch$family))==1,
                                  unique(wordMatch$family),
                                  NA)
      taxRow[1,"genus"] = ifelse(length(unique(wordMatch$genus))==1,
                                 unique(wordMatch$genus),
                                 NA)
      taxRow[1,"species"] = ifelse(length(unique(wordMatch$specificEpithet))==1,
                                   unique(wordMatch$specificEpithet),
                                   NA)
      return(taxRow)
    }
  }
  return(taxRow)
  
}))

taxTable = cbind(taxa.table, taxMatch)
write.csv(taxTable, "./rawdata/WFOTaxaTable.csv")
# some taxa require manual adjustment, which has been performed by authors
# on the EDITED version of this file.

# read in corrected taxa
taxTable = read.csv("./rawdata/WFOTaxaTableEDITED.csv")

# now apply corrected taxa to raw data

neotomaSyn = merge(plantSamples, taxTable[,c("taxa", "family", "genus", "species")],
                   by.x="variablename", by.y="taxa", all.x=TRUE, all.y=FALSE,
                   sort=FALSE)

neotomaSyn = neotomaSyn[order(neotomaSyn$siteid, neotomaSyn$age, neotomaSyn$family),]

neotomaSyn$family[is.na(neotomaSyn$family) | neotomaSyn$family == ""] = NA
neotomaSyn$genus[is.na(neotomaSyn$genus) | neotomaSyn$genus == ""] = NA
neotomaSyn$species[is.na(neotomaSyn$species) | neotomaSyn$species == ""] = NA

saveRDS(neotomaSyn, "./rawdata/processedRecords.rds")
