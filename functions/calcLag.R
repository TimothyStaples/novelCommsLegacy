calcLag <- function(novel.list, 
                            env.data, 
                            env.var, 
                            local.env.data,
                            local.lag,
                            global.lag){
  
# Prepare modelling data ####
print("Prepping data")

# aggregate 50 year temp change estimates to 200 year sampling bins, take mean
colnames(env.data)[colnames(env.data) == env.var] = "temp"

# global data
env.agg <- bin.env.data(env.data,
                        env.var = "temp",
                        bin.width=200,
                        lims=c(-100,25000))
env.agg <- env.agg[order(env.agg$bin, decreasing=TRUE), ]
  
# now difference temps based on required lag
env.lag <- diff(env.agg$env, differences=1, lag=unlist(global.lag))

env.agg$diff.env <- c(rep(NA, global.lag), env.lag)

local.lag.df <- t(sapply(split(local.env.data, f=local.env.data$siteid), function(x){
  x = x[order(as.numeric(x$bin), decreasing=TRUE),]
  c(rep(NA, local.lag), diff(x$localTemp, differences=1, lag=unlist(local.lag)))
}))
dimnames(local.lag.df) = list(levels(as.factor(local.env.data$siteid)), sort(as.numeric(unique(local.env.data$bin)), decreasing=TRUE))

# add in local temp data
local.lag.long <- long_form(dataTable = local.lag.df,
                            data.cols = matrix(rownames(local.lag.df), ncol=1, dimnames=list(rownames(local.lag.df), "siteid")),
                            category.cols = local.lag.df)
colnames(local.lag.long) = c("siteid", "bin", "local.diff")

comb.lag <- merge(local.lag.long, env.agg[,c("bin","diff.env")],
                  by.x="bin", by.y="bin",
                  all.x=TRUE, all.y=FALSE, sort=FALSE)
colnames(comb.lag) = c("bin", "siteid", paste0("localLag", local.lag),  paste0("globalLag",global.lag))

return(comb.lag)

}