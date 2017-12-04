originmodel <- function(variable_number){
  tmp.text1 = "V" %>% 
    paste(1:variable_number,sep = "") %>% 
    paste(",        lambda",sep = "") %>% 
    paste(1:variable_number,sep="") %>%
    paste(",    NA\n",sep = "")
  tmp.text1 = paste("LATFAC1 --> ",tmp.text1,sep="")
  tmp.text1 = paste(tmp.text1,collapse = "")
  tmp.text1 = paste("LATFAC1 <-> LATFAC1,   NA,          1\n", tmp.text1 ,sep="",collapse = "")
  tmp.text2 = "V" %>%
    paste(1:variable_number,sep = "") %>%
    paste(" <-> ",sep = "") %>%
    paste("V",sep = "") %>%
    paste(1:variable_number,sep = "") %>%
    paste(",             theta",sep = "") %>%
    paste(1:variable_number,sep = "") %>%
    paste(",     NA\n", sep = "") %>%
    paste(collapse = "") 
  tmp.text = paste(tmp.text1,tmp.text2,sep="",collapse = "")
  return(tmp.text)
}

altermodel <- function(origin_structure,whichvariable){
  tmp.text <- sub("NA,          1","phi11,       NA",origin_structure)
  tmp.text.alt <- sub(paste(c("lambda",whichvariable),collapse=""),"xyz",tmp.text) 
  tmp.text.alt <- sub("xyz,    NA","NA,          1",tmp.text.alt)
  return(tmp.text.alt)
}

cr.alpha = function(mmat){
  K = dim(mmat)[[2]]
  sigx = var(rowSums(mmat))
  sigy = sum(apply(mmat,2,var))
  return((K/(K-1))*(1-sigy/sigx))
}

GenerateData = function(n.sample,n.man,lambdas,gammas,taus){
  xi = scale(matrix(rnorm(n.sample),nrow = n.sample))
  noi = scale(matrix(rnorm(n.sample),nrow = n.sample))
  
  delta = scale(matrix(rnorm(n.sample*n.man),nrow = n.sample,ncol = n.man))
  
#  man.mean <- runif(n.man,min=2,max=6)
#  man.vars <- 1.6 - 0.5*abs(man.mean-4)
  XX.mat = array(0, dim = c(n.sample, n.man))
  pre = xi %*% lambdas + noi %*% gammas + delta %*% diag(sqrt(1 - lambdas^2 - gammas^2)) #standard matrix
  
#  pre = pre %*% diag(sqrt(man.vars)) + matrix(1,nrow=n.sample,ncol=1) %*% matrix(man.mean,nrow=1,ncol=n.man) #shifted matrix
#  pre = round(pre)
#  pre[which(pre>7)] = 7
#  pre[which(pre<1)] = 1 # boundary
  for(i in 1:length(taus)){ XX.mat = XX.mat + {pre < taus[i]} * 1 }
  return(XX.mat)
}

ContaminateData = function(DataMatrix, conta.percent){
  m = dim(DataMatrix)[2]
  n = dim(DataMatrix)[1]
  for(i in 1:round(conta.percent * n)){
    tmprow = sample(1:7,1)
    for(j in 2:m){
      tmprow = c(tmprow,sample(c(1:7)[-tmprow[j-1]],1))
    }
    DataMatrix[i,] = tmprow
  }
  return(DataMatrix)
}

generate.seed = function(n){
  my.new.seeds <- NULL
  for(i in 1:n) {
    tmp.ind <- floor(runif(1,min=2,max=627))
    next.rand <- .Random.seed[tmp.ind]
    my.new.seeds <- c(my.new.seeds,next.rand)
  }
  return(my.new.seeds)
}

load("names.rdata")

TestOriginModel = function(ctr, DataMatrix,semmodel,model.text){
  converge.tag = -1
  m = dim(DataMatrix)[2]
  DataFrame = as.data.frame(DataMatrix)
  ##################################################################  ################################################################## Step1
  do.chk = tryCatch({
    model.output = sem(semmodel,cov(DataFrame),dim(DataFrame)[1],debug = F, maxiter = 10000)
    summary.output = summary(model.output, fit.indices=c("RMSEA","samp.nFI","NNFI","CFI","SRMR"))
    do.chk = FALSE
  },error = function(e){TRUE})
  if(!do.chk){do.chk = (model.output$criterion < 0) || (model.output$iter > 9999) || (length(which(diag(model.output$vcov)<0)) > 0)}
  ##################################################################  ################################################################## Step2
  if(!do.chk){converge.tag = "Y"}
  else {
    nextModel = 1
    while(nextModel < m && do.chk){
      alter.text = altermodel(model.text, nextModel)
      model.sem.alter = specifyModel(text = alter.text,quiet = T)
      ##################################################################  ################################################################## 
      alt.chk = tryCatch({
        alter.output = sem(model.sem.alter, cov(DataFrame),dim(DataFrame)[1],debug = F, maxiter = 10000)
        summary.alter.output = summary(alter.output, fit.indices=c("RMSEA","samp.nFI","NNFI","CFI","SRMR"))
        alt.chk = FALSE
      },error = function(e){TRUE})
      if(!alt.chk){alt.chk = (alter.output$criterion < 0) || (alter.output$iter > 9999) || (length(which(diag(alter.output$vcov)<0)) > 0)}
      ##################################################################  ##################################################################
      if(!alt.chk){
        converge.tag = paste("ANCHOR",nextModel,sep = "")
        model.sem.tmp = semmodel
        model.sem.tmp[,3][-1][-nextModel] = alter.output$coeff[-1]
        ##################################################################  ##################################################################
        do.chk = tryCatch({
          tmp.output = sem(model.sem.tmp, cov(DataFrame),dim(DataFrame)[1],debug = F, maxiter = 10000)
          summary.tmp = summary(tmp.output, fit.indices=c("RMSEA","samp.nFI","NNFI","CFI","SRMR"))
          do.chk = FALSE
        },error = function(e){TRUE})
        if(!do.chk){do.chk = (tmp.output$criterion < 0) || (tmp.output$iter > 9999) || (length(which(diag(tmp.output$vcov)<0)) > 0)}
        ##################################################################  ##################################################################
        if(!do.chk){
          converge.tag = paste("YAn",nextModel,sep = "")
          eval(parse(text = paste("write.csv(x = DataMatrix, file = ","\"bad/orig_",ctr,".csv\")",sep = "")))
          model.output = tmp.output
          summary.output = summary.tmp
        }
      }
      nextModel = nextModel + 1
    }
  }
  if(!do.chk){ # get output if good model
    std.output = stdCoef(model.output)
    tmprow = data.frame(converge.type = converge.tag,
                        F.min.orig = model.output$criterion,
                        Chisq.orig = summary.output$chisq,
                        ChisqNull.orig = summary.output$chisqNull,
                        CFI.orig = summary.output$CFI,
                        NNFI.orig = summary.output$NNFI,
                        RMSEA.orig = summary.output$RMSEA[1],
                        SRMR.orig = summary.output$SRMR,
                        Df.orig = summary.output$df,
                        DfNULL.orig = summary.output$dfNull,
                        UCFI.orig = {(summary.output$chisqNull - summary.output$dfNull)-(summary.output$chisq - summary.output$df)}/(summary.output$chisqNull - summary.output$dfNull)
    )
    modelcoef = as.data.frame(t(c(std.output$`Std. Estimate`[-1][1:m], rep(0,12-m), std.output$`Std. Estimate`[-1][(1 + m):(2*m)], rep(0,12-m))))
    names(modelcoef) = c(paste("lambda",1:12,sep = ""),paste("theta",1:12,sep = ""))
    tmprow = cbind(tmprow, modelcoef)
    return(tmprow)
  }
  else{
    converge.tag = "N"
    tmprow = data.frame(converge.type = converge.tag)
    tmprow = cbind(tmprow, as.data.frame(t(rep(NA,34))))
    names(tmprow) = var.names1
    eval(parse(text = paste("write.csv(x = DataMatrix, file = ","\"bad/orig_",ctr,".csv\")",sep = "")))
    return(tmprow)
  }
}

TestContainmationModel = function(ctr, DataMatrix,semmodel,model.text){
  contamination.converge.tag = -1
  m = dim(DataMatrix)[2]
  DataFrame = as.data.frame(DataMatrix)
  ##################################################################  ################################################################## Step1
  do.chk = tryCatch({
    model.output = sem(semmodel,cov(DataFrame),dim(DataFrame)[1],debug = F, maxiter = 10000)
    summary.output = summary(model.output, fit.indices=c("RMSEA","samp.nFI","NNFI","CFI","SRMR"))
    do.chk = FALSE
  },error = function(e){TRUE})
  if(!do.chk){do.chk = (model.output$criterion < 0) || (model.output$iter > 9999) || (length(which(diag(model.output$vcov)<0)) > 0)}
  ##################################################################  ################################################################## Step2
  if(!do.chk){contamination.converge.tag = "Y"}
  else {
    nextModel = 1
    while(nextModel < m && do.chk){
      alter.text = altermodel(model.text, nextModel)
      model.sem.alter = specifyModel(text = alter.text,quiet = T)
      ##################################################################  ################################################################## 
      alt.chk = tryCatch({
        alter.output = sem(model.sem.alter, cov(DataFrame),dim(DataFrame)[1],debug = F, maxiter = 10000)
        summary.alter.output = summary(alter.output, fit.indices=c("RMSEA","samp.nFI","NNFI","CFI","SRMR"))
        alt.chk = FALSE
      },error = function(e){TRUE})
      if(!alt.chk){alt.chk = (alter.output$criterion < 0) || (alter.output$iter > 9999) || (length(which(diag(alter.output$vcov)<0)) > 0)}
      ##################################################################  ##################################################################
      if(!alt.chk){
        contamination.converge.tag = paste("ANCHOR",nextModel,sep = "")
        model.sem.tmp = semmodel
        model.sem.tmp[,3][-1][-nextModel] = alter.output$coeff[-1]
        ##################################################################  ##################################################################
        do.chk = tryCatch({
          tmp.output = sem(model.sem.tmp, cov(DataFrame),dim(DataFrame)[1],debug = F, maxiter = 10000)
          summary.tmp = summary(tmp.output, fit.indices=c("RMSEA","samp.nFI","NNFI","CFI","SRMR"))
          do.chk = FALSE
        },error = function(e){TRUE})
        if(!do.chk){do.chk = (tmp.output$criterion < 0) || (tmp.output$iter > 9999) || (length(which(diag(tmp.output$vcov)<0)) > 0)}
        ##################################################################  ##################################################################
        if(!do.chk){
          contamination.converge.tag = paste("YAn",nextModel,sep = "")
          eval(parse(text = paste("write.csv(x = DataMatrix, file = ","\"bad/",ctr,".csv\")",sep = "")))
          model.output = tmp.output
          summary.output = summary.tmp
        }
      }
      nextModel = nextModel + 1
    }
  }
  if(!do.chk){ # get output if good model
    std.output = stdCoef(model.output)
    tmprow = data.frame(conta.converge.type = contamination.converge.tag,
                        F.min = model.output$criterion,
                        Chisq = summary.output$chisq,
                        ChisqNull = summary.output$chisqNull,
                        CFI = summary.output$CFI,
                        NNFI = summary.output$NNFI,
                        RMSEA = summary.output$RMSEA[1],
                        SRMR = summary.output$SRMR,
                        Df = summary.output$df,
                        DfNULL = summary.output$dfNull,
                        UCFI = {(summary.output$chisqNull - summary.output$dfNull)-(summary.output$chisq - summary.output$df)}/(summary.output$chisqNull - summary.output$dfNull)
    )
    modelcoef = as.data.frame(t(c(std.output$`Std. Estimate`[-1][1:m], rep(0,12-m), std.output$`Std. Estimate`[-1][(1 + m):(2*m)], rep(0,12-m))))
    names(modelcoef) = c(paste("lambda",1:12,sep = ""),paste("theta",1:12,sep = ""))
    tmprow = cbind(tmprow, modelcoef)
    return(tmprow)
  }
  else{
    contamination.converge.tag = "N"
    tmprow = data.frame(conta.converge.type = contamination.converge.tag)
    tmprow = cbind(tmprow, as.data.frame(t(rep(NA,34))))
    names(tmprow) = var.names
    eval(parse(text = paste("write.csv(x = DataMatrix, file = ","\"bad/",ctr,".csv\")",sep = "")))
    return(tmprow)
  }
}


