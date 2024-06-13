create.ssmat <- function(records,
site.name, bins, age.limits = NA, tax.res, 
type = "pa" # options are "pa" (presence absence),
            #             "relpa (presence absence relativized by sampling effort)
            #             "abund" (relative abundance on subset of sites with abundance info)
){
  
  site.df <- droplevels(records[records$site == site.name &
                                !is.na(records$age), ])
  
  if(dim(site.df)[1]==0){return(NULL)}
  
  if(!is.na(age.limits[1])){
    
    site.df <- site.df[site.df$age <= age.limits[2] &
                         site.df$age >= age.limits[1], ]
    
  }
  
  # for each bin width
  site.df$age.bin <- cut(site.df$age, breaks = bins)

  age.bins <- levels(site.df$age.bin)
  mid.points <- cbind(age.bins,
                      as.numeric(substr(age.bins, regexpr(",", age.bins)+1,
                                        nchar(age.bins)-1)) - 0.5*(bins[2]-bins[1]))
  
  site.df$mid.point <- mid.points[match(site.df$age.bin, age.bins), 2]

  # cut out samples without taxonomic IDs
  site.sp <- droplevels(site.df[!is.na(site.df[, tax.res]), ])
  
  # filter out single observation sites
  if(length(unique(site.sp$age.bin[!is.na(site.sp$age.bin)]))<2 |
     length(unique(site.sp[!is.na(site.sp[,tax.res]),tax.res]))<2){
    return(NULL)
    }
  
  # cut out samples without species IDs & make site-species matrix
  if(tax.res == "species"){
    site.sp$binom <- paste(site.sp$genus, site.sp$species)
    site.sp.mat <- as.matrix(table(site.sp$age.bin, site.sp$binom))
    rownames(site.sp.mat) <- mid.points[match(rownames(site.sp.mat),
                                              age.bins), 2]
    } else {
    
    site.sp.mat <- as.matrix(table(site.sp$age.bin, site.sp[,tax.res]))
    site.sp.mat <- site.sp.mat[,colSums(site.sp.mat)>0]
    rownames(site.sp.mat) <- mid.points[match(rownames(site.sp.mat),
                                              age.bins), 2]
    
    }
  
  if(type=="pa"){
    
    site.sp.mat <- ifelse(site.sp.mat > 0, 1, 0)
    
  }
  
  if(type=="relpa"){
    
    temp <- site.df$age[site.df$mid.point== site.df$mid.point[1]]

    age.bin.rep <- tapply(site.df$age, site.df$mid.point, 
                          function(x){length(unique(x))})
    
    if(sum(names(age.bin.rep) %in% rownames(site.sp.mat)) != dim(site.sp.mat)[1]){
      print(paste0("ERROR in relativizing PA data in site",
                   site.df$site[1]))
      return(NULL)}
    
    site.sp.mat <- site.sp.mat / as.vector(age.bin.rep[match(rownames(site.sp.mat),
                                                             names(age.bin.rep))])
    
  }
  
  if(type == "abund"){
    
    site.sp$abund <- as.numeric(as.character(site.sp$count))
    
    site.sp$age.mid <- mid.points[match(site.sp$age.bin,
                                           age.bins),2]
    
    site.sp.mat <- matrix(0, 
                          nrow=sum(!is.na(unique(as.numeric(site.sp$age.mid)))),
                          ncol=length(unique(site.sp[,tax.res])))
    rownames(site.sp.mat) <- sort(unique(as.numeric(site.sp$age.mid)))
    colnames(site.sp.mat) <- unique(site.sp[,tax.res])
    
    for(tax in unique(site.sp[,tax.res])){

    temp <- site.sp[site.sp[,tax.res]==tax,] 

    tax.sum <- tapply(temp$abund, temp$age.mid, sum)
    
    site.sp.mat[names(tax.sum),as.character(tax)] = tax.sum
      
    }
    
  }
    
  site.sp.mat <- site.sp.mat[, colSums(site.sp.mat) > 0]
  
  return(as.matrix(unclass(site.sp.mat)))
  
}

