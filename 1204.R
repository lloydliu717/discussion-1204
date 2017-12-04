library(microbenchmark)

# lapply
lapply(1:3, function(x) c(x, x^2, x^3))
lapply(1:3/3, round, digits=3)

# parallel applys
library(parallel)
(no_cores = detectCores() - 1)

#The parallel package is basically about doing the above in parallel. 
#The main difference is that we need to start with setting up a cluster, a collection of “workers” that will be doing the job. 
#A good number of clusters is the numbers of available cores – 1.
#I’ve found that using all 8 cores on my machine will prevent me from doing anything else 
#(the computers comes to a standstill until the R task has finished). Therefore,

cl <- makeCluster(no_cores)
#stopCluster(cl)

parLapply(cl, 2:4,
          function(exponent) 2^exponent)

parLapply(cl, 1:3, function(x) c(x, x^2, x^3))

# variable scope
base <- 2

# fail
parLapply(cl, 2:4, 
          function(exponent) base^exponent)

clusterExport(cl, "base")

# pass
parLapply(cl, 2:4, 
          function(exponent) base^exponent)
# once you broadcast the variable, that any changes to the variable after clusterExport are ignored
base = 3

parLapply(cl, 2:4, 
          function(exponent) base^exponent)

# sapply
parSapply(cl, 2:4, 
          function(exponent) 
            base^exponent)

parSapply(cl, 2:4, 
          function(exponent){
            c(base = base^exponent, self = exponent^exponent)
          })

stopCluster(cl)


# foreach
library(foreach)
library(doMC)
#library(doParallel)

registerDoMC(no_cores)
#cl<-makeCluster(no_cores)
#registerDoParallel(cl)


getDoParWorkers()

foreach(exponent = 2:4, .combine = c)  %dopar%  base^exponent # c/cbind/rbind/list

foreach(exponent = 2:4, .combine = list, .multicombine = TRUE)  %dopar%  base^exponent
#the .multicombine argument that is needed to avoid a nested list. 
#The nesting occurs due to the sequential .combine function calls, 
#i.e. list(list(result.1, result.2), result.3):

foreach(exponent = 2:4, .combine = , .export = )

# simulate some data
D <- rnorm(1000, 165, 5)
M <- D + rnorm(1000, 0, 1)

# calculate linear model
DataModel <- lm(D ~ M)
Beta <- coef( DataModel )[2]

nSim<-10000

# Execute sampling and analysis in parallel
matrix <- foreach(i=1:nSim, .combine=rbind) %dopar% {
  perm <- sample(D, replace=FALSE)
  mdl <- lm(perm ~ M)
  c(i, coef(mdl))
}

