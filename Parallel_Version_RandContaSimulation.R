library(doMC)           # needed for parallel computing
library(sem)            # needed for latent variable analysis
library(parallel)       # needed for parallel computing
library(foreach)        # needed for parallel computing
library(iterators)      # needed for parallel computing
#setwd("/Users/liuzichun/Documents/BU/2017SPRING/Characteristic_Respondents/SurveyDataSimulation/RandContaSimulation/")

#####detect cores and register multiple cores to work
#num_cores <- detectCores()       ## Calculate the number of available cores

num_cores = as.numeric(Sys.getenv("NSLOTS"))

if(is.na(num_cores)){num_cores = 1}

# num_cores <- num_cores - 1       ## always best to use one less than max when using MAC
max_cores <-  28                  ## set to a given value or use 8 or 16
##     on the cluster this value whould be set to the
##     number of cores requested by the job
##     on scc-lite this ncores should be set to any number
##     between 1 and 28 or 36
# max_cores <- 3

num_cores <- min(max_cores,num_cores)

registerDoMC(num_cores)

getDoParWorkers()



######some functions
source("Function1.R")
# model text generation
# alternative model text generation
# calculate cr.alpha
# generate dataset
# generate corresponding contaminated dataset
# test model on original dataset and error processing
# test model on contaminated dataset and error processing
# generate a series of random seed

######set random seed
#my.new.seeds = generate.seed(25)

###### global value setup
source("Setup.R")
# granularity 
# minimum lambda
# number of manifest items
# sample size
# percentage contamination
# pure repeat times

###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
###### ###### ###### ###### ######
###### ###### ###### ###### ######
###### ###### ###### ###### ###### Parallel loop
ptm = proc.time()
ctr = 0
MDF = data.frame()

print(Sys.time())
print("Start!!!!!!!!")
MDF = foreach(l = 1:n.runs,.combine = rbind) %dopar% {
  for(i in 1:length(n.man)){
    ori.model.text = originmodel(n.man[i])
    model.sem.SIM = specifyModel(text = ori.model.text,quiet = T)
    for(g in gammalist){
      gammas = g * {{-1}^{1:n.man[i]}+1}/2
      for(minla in minlambda){
        lambdas = seq(0.86, minla, length.out = n.man[i])
        for(ta in taubound){
          taus = seq(-ta,ta,length.out = 6)
          for(j in 1:length(n.sample)){
            for(k in 1:length(per.conta)){
              ctr = ctr + 1
              #"Generate original data set"
              Xmat.orig = GenerateData(n.sample = n.sample[j], n.man = n.man[i],lambdas,gammas,taus)
              #"Generate conta data set"
              Xmat.conta = ContaminateData(DataMatrix = Xmat.orig, conta.percent = per.conta[k])
              #"Record loop and data set information"
              container1 = data.frame(ctr = ctr, 
                                      manifest.item = n.man[i], 
                                      sample.size = n.sample[j], 
                                      min.lambda = minla,
                                      gamma = g,
                                      taubound = ta, 
                                      noise.range = paste("0","0.4",sep = ","),
                                      contamination.percent = per.conta[k],
                                      cr.alpha.orig = cr.alpha(Xmat.orig), 
                                      cr.alpha.conta = cr.alpha(Xmat.conta))
              #"run original sem model"
              container2 = TestOriginModel(ctr,Xmat.orig, model.sem.SIM, ori.model.text)
              #"run conta sem model"
              container3 = TestContainmationModel(ctr, Xmat.conta, model.sem.SIM, ori.model.text)
              
              MDF = rbind(MDF, cbind(container1, container2, container3))
            }
          }
        }
      }
    }
  }  
  return(MDF)
}
proc.time() - ptm
print(paste("In total get",nrow(MDF),"rows"))
write.csv(MDF,"MDF_4_14.csv")
