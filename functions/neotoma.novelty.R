neotoma.novelty <- function(dataset,
                            ssmat.type = "pa",
                            taxon.res,
                            bins,
                            age.limits,
                            bin.cutoff,
                            taxa.cutoff,
                            novel.alpha,
                            novel.metric,
                            rich.cutoff,
                            novel.plot = TRUE,
                            sqrt.mat=FALSE){
  
  print("Creating site-species matrices...")
  # create site-species matrices using bin width

  cand.ssmats <- lapply(unique(dataset$site), function(site){
    print(site)
    temp <- create.ssmat(records = dataset,
                 site.name = site,
                 bins = bins,
                 tax.res = taxon.res,
                 type=ssmat.type,
                 age.limits = age.limits)
    
    if(class(temp)[1]=="matrix"){
  
    # remove bins within sites with fewer counts than cut off
    temp <- temp[rowSums(temp) >= rich.cutoff[1] &
                 rowSums(temp) <= rich.cutoff[2],]

    }
    
    return(temp)
    
  })
  names(cand.ssmats) <- unique(dataset$site)
  
  # remove sites that had no appropriate data
  cand.ssmats <- cand.ssmats[sapply(cand.ssmats, function(x){class(x)[1]})=="matrix"]
  
  # remove sites with fewer age bins or taxa than cutoff thresholds
  ssmat.dims <- t(sapply(cand.ssmats, dim))
  cand.ssmats <- cand.ssmats[ssmat.dims[,1] >= bin.cutoff & 
                             ssmat.dims[,2] >= taxa.cutoff]
  
  # do we square-root transform first?
  if(sqrt.mat){
  cand.ssmats <- lapply(cand.ssmats, sqrt)
  }

  cand.prop <- lapply(cand.ssmats, function(x){prop.table(x, 1)})
  
  # use the identify.novel.gam function to find novel communities in each time-series.
  # function creates plots as it runs.
  print("Identifying novel communities...")
  cand.novel <- lapply(names(cand.prop), function(site){
      identify.novel.gam(site.sp.mat = cand.prop[[site]], 
                         alpha = novel.alpha, 
                         metric=novel.metric,
                         site=site,
                         plot = novel.plot)
  })
  names(cand.novel) <- names(cand.ssmats)
  
  # remove communities that couldn't be estimated due to no variation
  # in taxa
  cand.subset <- sapply(cand.novel, class) == "data.frame"
  cand.ssmats <- cand.ssmats[cand.subset]
  cand.novel <- cand.novel[cand.subset]
  
  return(list(data = dataset[dataset$site %in% names(cand.ssmats), ],
              sites = names(cand.ssmats),
              raw.ssmats = cand.ssmats,
              prop.ssmats = cand.prop,
              novel = cand.novel))
  
}
