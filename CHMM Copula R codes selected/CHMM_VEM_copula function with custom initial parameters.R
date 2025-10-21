res.init.esProb<-matrix(c(0.3,0.4,0.45,0.5,0.55,0.65),byrow = T,2,3)
res.init.transPr<-p3
res.init.initPr<-c(0.1,0.8,0.1)

###-------------------------------------------------
CHMM_VEM <- function(X, nb.states, 
                     itmax = 5e2,
                     threshold = 1e-7,nbI,nbD,omega,res.init.esProb,res.init.transPr,res.init.initPr){ 
  nbT <- nrow(X) 
  
  
  
  #initial parameters
  res.init <- init.VEM(X=X, nb.states=nb.states, nbD=nbD, nbT=nbT, nbI=nbI)
  esProb <- res.init.esProb  
  transPr <- res.init.transPr
  initPr <- res.init.initPr
  postPr.last <- res.init$postPr
  
  #
  Old.param = c(esProb, transPr)
  for (iter in 1:itmax) {
    postPr.list <- list()
    emisPr.list <- list()
    emisPrW.list <- list()
    trsTmp.list<-list()
    trsTmp <- matrix(0, nb.states, nb.states)
    emisPr <- matrix(0, nbT, nb.states)
    
    
    #---------------------------------
    copula<-matrix(0,nbD,nbD)
    sumlogcopula<-rep(NA,nbD)
    for (dis in 1:nbD) {
      for (dis2 in c(1:nbD)[-dis]){ 
        om<-omega[dis,dis2]
        binomial_copula=matrix(NA,nb.states,nb.states) #n is the number of trials
        binomial_copula[1,1]<- 1/3*om*(om+1)/(om^2+om+1+sqrt(om*(om+2)*(2*om+1)))
        binomial_copula[1,2]<- 1/3*(sqrt(om*(om+2)*(2*om+1))-3*om)/((om-1)^2)
        binomial_copula[1,3]<- 1/3*(om+1)/(om^2+om+1+sqrt(om*(om+2)*(2*om+1)))
        binomial_copula[2,1]<- 1/3*(sqrt(om*(om+2)*(2*om+1))-3*om)/((om-1)^2)
        binomial_copula[2,2]<- 1/3*(om^2+4*om+1-2*sqrt(om*(om+2)*(2*om+1)))/((om-1)^2)
        binomial_copula[2,3]<- 1/3*(sqrt(om*(om+2)*(2*om+1))-3*om)/((om-1)^2)
        binomial_copula[3,1]<- 1/3*(om+1)/(om^2+om+1+sqrt(om*(om+2)*(2*om+1)))
        binomial_copula[3,2]<- 1/3*(sqrt(om*(om+2)*(2*om+1))-3*om)/((om-1)^2)
        binomial_copula[3,3]<- 1/3*om*(om+1)/(om^2+om+1+sqrt(om*(om+2)*(2*om+1)))
        
        copula[dis,dis2]<- max(binomial_copula)
        
      }
      
      sumlogcopula[dis]<-sum(log(copula[dis,-dis]))  
      
      
      for (t in 1:nbT)  {
        trsTmp <- matrix(0, nb.states, nb.states)
        
        emisPr[t,] <- dbinom(X[t,dis],size=nbI[dis],prob=esProb[dis,]) #for each disease
        emisPr.list[[dis]] <- emisPr 
        emisPrW <- emisPr*sumlogcopula[dis] 
        emisPrW.list[[dis]] <- emisPrW           #here we have h_it for each disease and state
        
        
        # Transform matrix to vectors     -----------------------------
        initTmp <- initPr * emisPrW[1,]  
        initPr <- initTmp /sum(initTmp)             #normalized
        
        # Forward-Backward recursion      -----------------------------
        resF <- ForwardR(initPr, emisPrW, transPr, nbT, nb.states)
        resB <- BackwardR(resF$Fpr, transPr, nbT, nb.states)
        
        
        
        postPr.tmp <- resB$postPr
        #postPr.tmp <- apply(postPr.tmp, 2, pmax, 1e-200)
        postPr <- postPr.tmp / rowSums(postPr.tmp)
        postPr.list[[dis]] <- postPr
        Fpr <- resF$Fpr
        Gpr <- resB$Gpr
        
        trsTmp <- trsTmp+transPr * t(Fpr[-nbT,]) %*% (postPr[-1,] / Gpr[-1,])
        
      }
      trsTmp.list[[dis]]<-trsTmp
    }
    
    postPr.last <- postPr.list
    
    ## 2.M-step      ----------------------------
    #-- for aggregate individuals calculate transit, init, emis parameters
    # update transPr
    
    trsAvg<-Reduce("+",trsTmp.list)
    transPr <- trsAvg / rowSums(trsAvg)
    
    # update initPr       ----------------------------
    init.tmp <- rep(0, nb.states)
    for(i in 1:nbD) init.tmp <- init.tmp + postPr.list[[i]][1,]
    initPr <- init.tmp / sum(init.tmp)
    
    for (dis in 1:nbD){
      esProb[dis,]<-colSums(postPr.list[[dis]]*X[,dis])/nbI[dis]/colSums(postPr.list[[dis]])
    }
    
    
    #  order the parameters          ----------------------------
    #ordProb <- order(esProb)
    #esProb <- sort(esProb)
    #transPr <- transPr[ordProb,ordProb]
    
    ## 3.Stop iteration          ----------------------------
    New.param <-  c(esProb, transPr)
    crit <- New.param - Old.param
    Old.param  <- New.param
    
    if (iter > 1 && max(abs(crit)) <= threshold) break()
  }
  return(list(res.init=res.init,postPr = postPr.list, initPr = initPr,  transPr = transPr,  esProb = esProb,
              emisPr = emisPr.list, emisPrW = emisPrW.list, iterstop = iter))
}




Model<-CHMM_VEM(X=X, nb.states=3, itmax = 5e2, threshold = 1e-7, nbI=nbI,nbD=nbD, omega=omega,res.init.esProb=res.init.esProb,
         res.init.transPr=res.init.transPr,res.init.initPr=res.init.initPr)
