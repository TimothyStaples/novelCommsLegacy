mer.ci <- function(model, newdata, sims, ncores=1, re.form=NA){
  
  if(!is.na(re.form)){re.form = as.formula(re.form)}
    # bootstrap predictions 
    bb <- bootMer(model,
                  nsim=sims,
                  FUN=function(fit){predict(fit, newdata, re.form=re.form)},
                  parallel="multicore",
                  ncpus = ncores)
  
  # extract quantiles over simulations to get CIs
  cis <- t(apply(bb$t, 2, function(x){quantile(x, c(0.025, 0.975), na.rm=TRUE)}))
  
  return(data.frame(fit = bb$t0,
                    lower = cis[,1],
                    upper = cis[,2]))
  
}
