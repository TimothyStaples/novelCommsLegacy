bin.env.data <- function(env.data, bin.width, lims, env.var){

env.data <- env.data[order(env.data$age, decreasing=TRUE),]
# env.data <- env.data[complete.cases(env.data),]

env.data$env.agg <- cut(env.data$age, breaks=seq(lims[1],lims[2],bin.width))

env.agg.char <- as.character(env.data$env.agg)
env.data$bin <- as.numeric(substr(env.agg.char, regexpr(",", env.agg.char)+1,
                                  nchar(env.agg.char)-1)) - 0.5*bin.width


env.data$env.var <- env.data[,env.var]

env.data <- with(env.data,
                 data.frame(env = tapply(env.var, bin, mean, na.rm=TRUE),
                            env.sd = tapply(env.var, bin, sd, na.rm=TRUE),
                            env.se = tapply(env.var, bin, function(x){
                              sd(x) / sqrt(length(x))}),
                            bin = as.numeric(levels(as.factor(bin)))))

return(env.data)
}