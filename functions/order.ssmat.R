order.ssmat <- function(site.sp.mat){
  
  if(as.numeric(as.character(rownames(site.sp.mat)[1])) < 
     as.numeric(as.character(rownames(site.sp.mat)[dim(site.sp.mat)[1]]))){
    
    site.sp.mat <- site.sp.mat[dim(site.sp.mat)[1]:1, ]
    
  }

return(site.sp.mat)

}