novel.trajectoryV4 <- function(novel.object,
                               ssmat.type,
                               dissim.method,
                               genusNov,
                               novel.crop = c(-Inf, Inf)){

  print("Triangulating observed trajectory")

  require(hilldiv)
  require(vegan)
  
  sites <- novel.object$sites
  novel.list <- novel.object$novel
  ssmat.list <- novel.object[[ssmat.type]]

    novel.traj.list <- lapply(c("instant", "cumul", "novel"), function(novel.cat){
    
    print(novel.cat)
      
    # use s = "3502" if you need an example for this lapply
    do.call("rbind", lapply(sites, function(s){
      
      print(s)
      temp.nov <- novel.list[[s]]
      temp.ssmat <- ssmat.list[[s]]
      temp.full = novel.object[["raw.ssmats"]][[s]]
      
      # find position of novel communities
      nov.locs <- as.character(temp.nov$bins)[temp.nov[,novel.cat]]
      
      # subset novel communities to just those that occur in certain time band
      nov.locs <- nov.locs[as.numeric(nov.locs) >= novel.crop[1] &
                           as.numeric(nov.locs) <= novel.crop[2]]
      nov.locs <- nov.locs[!is.na(nov.locs)]
      
      if(length(nov.locs) == 0){return(NULL)}
      
      # return data-frame of post-novel trajectories
      site.df <- do.call("rbind", lapply(nov.locs, function(l){
        
        print(paste0("--- ", l))
        # rotate ssmat if we need to
        if(as.numeric(as.character(rownames(temp.ssmat)[1])) < 
           as.numeric(as.character(rownames(temp.ssmat)[dim(temp.ssmat)[1]]))){
          
          temp.ssmat <- temp.ssmat[dim(temp.ssmat)[1]:1, ]
          temp.full <- temp.full[dim(temp.full)[1]:1,]
        }
        
        # only keep going if novel comm is at least 2 time points from time series end
        if(which(rownames(temp.ssmat) == l)+1 >= nrow(temp.ssmat)){return(NULL)}
        
        # subset ssmat from novel time point -1 ("pre-novel")
        sub.ssmat <- temp.ssmat[(which(rownames(temp.ssmat) == l)-1):nrow(temp.ssmat),]
        sub.dist <- as.matrix(vegdist(sub.ssmat, method=dissim.method))
        sub.raw = temp.full[(which(rownames(temp.full) == l)-1):nrow(temp.full),]
        
        # create data-frame of dissimilarities after novel comm
        dist.df <- data.frame(sub.dist[-(1:2),1:2])
        
        if(ncol(dist.df)==1){return(NULL)}
        
        colnames(dist.df) <- c("dP", "dN")
        dist.df$delta.dP <- diff(sub.dist[-1,1])
        dist.df$delta.dN <- c(NA, diff(sub.dist[-(1:2),2]))
        dist.df$dP1 <- dist.df$dP[1]
        dist.df$dN1 <- dist.df$dN[1]
        dist.df$bin <- as.numeric(rownames(dist.df))
        dist.df$novel.time <- as.numeric(l)
        dist.df$time.since.novel <- abs(dist.df$bin - as.numeric(l))
        dist.df$site <- s
        dist.df$novelID <- paste0(dist.df$site, ".", dist.df$novel.time)
        dist.df$bin.n <- which(rownames(temp.ssmat) == l)
        
        dist.df$H0 = apply(sub.ssmat[-(1:2),],1, function(x){hill_div(x[x>0], qvalue=0)})
        dist.df$H1 = apply(sub.ssmat[-(1:2),],1, function(x){hill_div(x[x>0], qvalue=1)})
        dist.df$H2 = apply(sub.ssmat[-(1:2),],1, function(x){hill_div(x[x>0], qvalue=2)})
        
        rawRound = t(apply(sub.raw[-(1:2),]*10, 1, round))
        dist.df$rareDivN = ifelse(min(rowSums(rawRound)) < 500, min(rowSums(rawRound)), 500)
        dist.df$rareType = rarefy(rawRound, 10)
        dist.df$rareDiv = rarefy(rawRound, dist.df$rareDivN[1])
        
        novelComm <- sub.ssmat[2,][sub.ssmat[2,] > 0]
        dist.df$novelH0 = hill_div(novelComm, qvalue=0)
        dist.df$novelH1 = hill_div(novelComm, qvalue=1)
        dist.df$novelH2 = hill_div(novelComm, qvalue=2)
        dist.df$novelRareType = rarefy(round(sub.raw[2,]*10), 10)
        dist.df$novelRareDiv = rarefy(round(sub.raw[2,]*10), dist.df$rareDivN[1])
        
        preComm <- sub.ssmat[1,][sub.ssmat[1,] > 0]
        dist.df$preNovH0 = hill_div(preComm, qvalue=0)
        dist.df$preNovH1 = hill_div(preComm, qvalue=1)
        dist.df$preNovH2 = hill_div(preComm, qvalue=2)
        dist.df$preRareType = rarefy(round(sub.raw[1,]*10), 10)
        dist.df$preRareDiv = rarefy(round(sub.raw[1,]*10), dist.df$rareDivN[1])
        
        # get pollen count of the bin (abundance effect)
        rawMat <- novel.object$raw.ssmats[[s]]
        dist.df$novAbund <- sum(rawMat[l,])
        dist.df$novAbundProp <- sum(rawMat[l,], na.rm=TRUE) / mean(rowSums(rawMat[rownames(rawMat) != l,], na.rm=TRUE),
                                                       na.rm=TRUE)
        
        # tease out increases and decreases
        tax.df = do.call('rbind', lapply(3:nrow(sub.ssmat), function(n){
          rowDiffs = sub.ssmat[n-1,] - sub.ssmat[n,]
          data.frame(downerN = sum(rowDiffs < 0), upperN = sum(rowDiffs > 0), relAbundChange = sum(rowDiffs[rowDiffs > 0]))
        }))
        
        dist.df = cbind(dist.df, tax.df)
        
        return(dist.df)
      }))
      
      if(is.null(site.df)){return(NULL)}
      
      site.df$gamma <- ncol(novel.object$raw.ssmats[[s]])
      site.df$tslength <- nrow(novel.object$raw.ssmats[[s]])
      
      return(site.df)
    }))
  
    })
    names(novel.traj.list) <- c("instant", "cumul", "novel")
    
    # now create the same object using randomly selected, non-novel comms
    print("Triangulated non-novel trajectory")
    
    rand.traj.df <- do.call("rbind", lapply(sites, function(s){
      #print(s)
      temp.nov <- novel.list[[s]]
      temp.ssmat <- ssmat.list[[s]]
      temp.full = novel.object[["raw.ssmats"]][[s]]
      
      # find position of novel communities
      back.locs <- as.character(temp.nov$bins)[temp.nov$cat == "back"]
      
      # subset novel communities to just those that occur in certain time band
      back.locs <- back.locs[as.numeric(back.locs) >= novel.crop[1] &
                             as.numeric(back.locs) <= novel.crop[2]]
      back.locs <- back.locs[!is.na(back.locs)]
      
      if(length(back.locs) == 0){return(NULL)}
      
      # return data-frame of post-novel trajectories  # return data-frame of post-novel trajectories
      site.df <- do.call("rbind", lapply(back.locs, function(l){
        
        print(paste0("--- ", l))
        # rotate ssmat if we need to
        if(as.numeric(as.character(rownames(temp.ssmat)[1])) < 
           as.numeric(as.character(rownames(temp.ssmat)[dim(temp.ssmat)[1]]))){
          
          temp.ssmat <- temp.ssmat[dim(temp.ssmat)[1]:1, ]
          temp.full <- temp.full[dim(temp.full)[1]:1,]
        }
        
        # only keep going if novel comm is at least 2 time points from time series end
        if(which(rownames(temp.ssmat) == l)+1 >= nrow(temp.ssmat)){return(NULL)}
        
        # subset ssmat from novel time point -1 ("pre-novel")
        sub.ssmat <- temp.ssmat[(which(rownames(temp.ssmat) == l)-1):nrow(temp.ssmat),]
        sub.dist <- as.matrix(vegdist(sub.ssmat, method=dissim.method))
        sub.raw = temp.full[(which(rownames(temp.full) == l)-1):nrow(temp.full),]
        
        # create data-frame of dissimilarities after novel comm
        dist.df <- data.frame(sub.dist[-(1:2),1:2])
        
        if(ncol(dist.df)==1){return(NULL)}
        
        colnames(dist.df) <- c("dP", "dN")
        dist.df$delta.dP <- diff(sub.dist[-1,1])
        dist.df$delta.dN <- c(NA, diff(sub.dist[-(1:2),2]))
        dist.df$dP1 <- dist.df$dP[1]
        dist.df$dN1 <- dist.df$dN[1]
        dist.df$bin <- as.numeric(rownames(dist.df))
        dist.df$novel.time <- as.numeric(l)
        dist.df$time.since.novel <- abs(dist.df$bin - as.numeric(l))
        dist.df$site <- s
        dist.df$novelID <- paste0(dist.df$site, ".", dist.df$novel.time)
        dist.df$bin.n <- which(rownames(temp.ssmat) == l)
        
        dist.df$H0 = apply(sub.ssmat[-(1:2),],1, function(x){hill_div(x[x>0], qvalue=0)})
        dist.df$H1 = apply(sub.ssmat[-(1:2),],1, function(x){hill_div(x[x>0], qvalue=1)})
        dist.df$H2 = apply(sub.ssmat[-(1:2),],1, function(x){hill_div(x[x>0], qvalue=2)})
        
        rawRound = t(apply(sub.raw[-(1:2),]*10, 1, round))
        dist.df$rareDivN = ifelse(min(rowSums(rawRound)) < 500, min(rowSums(rawRound)), 500)
        dist.df$rareType = rarefy(rawRound, 10)
        dist.df$rareDiv = rarefy(rawRound, dist.df$rareDivN[1])
        
        novelComm <- sub.ssmat[2,][sub.ssmat[2,] > 0]
        dist.df$novelH0 = hill_div(novelComm, qvalue=0)
        dist.df$novelH1 = hill_div(novelComm, qvalue=1)
        dist.df$novelH2 = hill_div(novelComm, qvalue=2)
        dist.df$novelRareType = rarefy(round(sub.raw[2,]*10), 10)
        dist.df$novelRareDiv = rarefy(round(sub.raw[2,]*10), dist.df$rareDivN[1])
        
        preComm <- sub.ssmat[1,][sub.ssmat[1,] > 0]
        dist.df$preNovH0 = hill_div(preComm, qvalue=0)
        dist.df$preNovH1 = hill_div(preComm, qvalue=1)
        dist.df$preNovH2 = hill_div(preComm, qvalue=2)
        dist.df$preRareType = rarefy(round(sub.raw[1,]*10), 10)
        dist.df$preRareDiv = rarefy(round(sub.raw[1,]*10), dist.df$rareDivN[1])
        
        # get pollen count of the bin (abundance effect)
        rawMat <- novel.object$raw.ssmats[[s]]
        dist.df$novAbund <- sum(rawMat[l,])
        dist.df$novAbundProp <- sum(rawMat[l,], na.rm=TRUE) / mean(rowSums(rawMat[rownames(rawMat) != l,], na.rm=TRUE),
                                                                   na.rm=TRUE)
        
        return(dist.df)
      }))
      
      if(is.null(site.df)){return(NULL)}
      
      site.df$gamma <- ncol(novel.object$raw.ssmats[[s]])
      site.df$tslength <- nrow(novel.object$raw.ssmats[[s]])
      
      return(site.df)
      
    }))
    
    return(list(novel.traj.list = novel.traj.list,
                rand.traj.list = rand.traj.df))  
    
}
