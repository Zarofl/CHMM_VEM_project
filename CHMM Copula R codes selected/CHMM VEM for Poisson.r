CHMM_VEM <- function(X, nb.states, 
                     itmax = 5e2,
                     threshold = 1e-7,nbD,omega){ 
  nbT <- nrow(X) 
  
  #initial parameters
  res.init <- init.VEM(X=X, nb.states=nb.states, nbD=nbD, nbT=nbT)
  esLambda <- res.init$esLambda  
  transPr <- res.init$transPr
  initPr <- res.init$initPr
  postPr.last <- res.init$postPr
  
  
  #
  Old.param = c(esLambda, transPr)
  for (iter in 1:itmax) {
    lpostPr.list <- list()
    emisPr.list <- list()
    emisPrW.list <- list()
    trsTmp.list<-list()
    resF.list<-list()
    resB.list<-list()
    llk.list<-list()
    lFpr.list<-list()
    lGpr.list<-list()
    trsTmp <- matrix(0, nb.states, nb.states)
    emisPr <- matrix(0, nbT, nb.states)
    AIC<-rep(NA,nbD)
    BIC<-rep(NA,nbD)
    
    
    #-----Calculate Binomial copula
    copula<-matrix(0,nbD,nbD)
    prodcopula<-rep(NA,nbD)
    RSS<-matrix(NA,nbD,nbD)
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
      
      prodcopula[dis]<-prod(copula[dis,-dis])  
      
      lpost_z<-rep(NA,nbT)
      
      trsTmp <- matrix(0, nb.states, nb.states) 
      
      #Calculate emission probabilitie that follow Poisson
      
      for (t in 1:nbT)  {
        emisPr[t,] <- dpois(X[t,dis],lambda=esLambda[dis,]) #for each disease 
      } 
      print(emisPr)
      
      emisPr.list[[dis]] <- emisPr 
      emisPrW <- emisPr*prodcopula[dis]       
      emisPrW.list[[dis]] <- emisPrW           
      
      
      # initPr are also affected by h_it
      initTmp <- initPr * emisPrW[1,]  
      initPr <- initTmp /sum(initTmp)  #normalized
      
      
      # Forward-Backward recursion      -----------------------------
      resF <- ForwardR(initPr, emisPrW, transPr, nbT, nb.states)
      resF.list[[dis]]<-resF$lFpr
      resB <- BackwardR(resF$lFpr, transPr, nbT, nb.states)

      #Here we have a problem with underflow, we change "-inf" to the mean value
      for (t in 1:nbT){
        for (i in 1:nb.states){
          if (resB$lpostPr[t,i]=="-Inf") {
            matnoinf<-as.data.frame(resB$lpostPr) %>% 
              filter(resB$lpostPr[,i]!="-Inf")
            resB$lpostPr[t,i]<-mean(matnoinf[,i])
          }
        }
      }
      
      resB.list[[dis]]<-resB
      
      
      lpostPr.tmp <- resB$lpostPr
      #lpostPr.tmp <- apply(lpostPr.tmp, 2, pmax, -100)     
      
      for (t in 1:nbT){
        
        lpost_z[t]<-max(lpostPr.tmp[t,])+log(sum(exp(lpostPr.tmp[t,])-max(lpostPr.tmp[t,]))) #can use mpfr here
      }  
      
      
      lpostPr <- lpostPr.tmp - lpost_z 
      
      
      # postPr <- postPr.tmp / rowSums(postPr.tmp)  #RSS are very high, model fit is low
      
      #lpostPr <- lpostPr.tmp / rowSums(lpostPr.tmp)
      
      lpostPr.list[[dis]] <- lpostPr
      lFpr <- resF$lFpr
      lGpr <- resB$lGpr
      
      #--Here we calculate matrix using formulas in backward recursion, actually we calculate conditional expectation
      # c<-asNumeric(exp(mpfr(lFpr[-nbT,],5)))
      c<-exp(mpfr(lFpr[-nbT,],5))
      
      #  d<-asNumeric(exp(mpfr(lpostPr[-1,] - lGpr[-1,],5)))
      d<-exp(mpfr(lpostPr[-1,] - lGpr[-1,],5))  
      
      trsTmp <- trsTmp+transPr * t(c) %*% d 
      trsTmp<-as.numeric(trsTmp)
      trsTmp<-matrix(trsTmp,nb.states,nb.states)
      
      #trsTmp <- trsTmp+transPr * t(exp(mpfr(lFpr[-nbT,],5))) %*% (exp(mpfr(lpostPr[-1,] - lGpr[-1,],5))) #do same with numeric
      
      
      trsTmp.list[[dis]]<-trsTmp
      lFpr.list[[dis]]<-lFpr
      lGpr.list[[dis]]<-lGpr
      
      
    }
    
    

    ## 2.M-step      ----------------------------
    #-- for aggregate individuals calculate transit, init, emis parameters
    # update transPr
    
    trsAvg<-Reduce("+",trsTmp.list)
    transPr <- trsAvg / rowSums(trsAvg)

    # update initPr       ----------------------------
    init.tmp <- rep(0, nb.states)
    
    for(i in 1:nbD){
      e<-as.numeric(exp(mpfr(lpostPr.list[[i]][1,],5)))
      e<-matrix(e,1,nb.states)
      
      init.tmp <- init.tmp + e
    }
    initPr <- init.tmp / sum(init.tmp)
    
    #update Lambdas for emission distribution----------------------
    
    for (dis in 1:nbD){
      c<-exp(mpfr(lFpr.list[[dis]][-nbT,],5))
      c<-asNumeric(c)

      d<-exp(mpfr(lpostPr.list[[dis]][-1,] - lGpr.list[[dis]][-1,],5))  
      d<-asNumeric(d)

      esLambda[dis,]<-colSums(X[-nbT,dis]*c*d)/colSums(c*d)
    }


    ## Calculate criteria    RSS     -------------------------- I will update it 
    
    # for (dis in 1:nbD) {
    #   for (dis2 in c(1:nbD)[-dis]){
    #     a<-as.numeric(exp(mpfr(lpostPr.list[[dis]],5)))
    #     a<-matrix(a,nbT,nb.states)
    #     b<-as.numeric(exp(mpfr(lpostPr.list[[dis2]],5)))
    #     b<-matrix(b,nbT,nb.states)
    #     RSS[dis,dis2]<-sum(a*(X[,dis]-nbI[dis]*esProb[dis,])^2)+sum(b*(X[,dis2]-nbI[dis2]*esProb[dis2,])^2)
    #     #RSS[dis,dis2]<-sum(exp(mpfr(lpostPr.list[[dis]],5))*(X[,dis]-nbI[dis]*esProb[dis,])^2)+sum(exp(mpfr(lpostPr.list[[dis2]],5))*(X[,dis2]-nbI[dis2]*esProb[dis2,])^2)
    #   }
    # }
    
    ## Calculate AIC and BIC
    
    for (dis in 1:nbD){
      c<-max(resF.list[[dis]][nbT,])
      llk.list[[dis]]<-c+log(sum(exp(resF.list[[dis]][nbT,]-c)))
    }
    
    numb.par<-nbD*(nbD-1)/2+nbD*nb.states+nb.states*(nb.states-1)+nb.states-1 
    
    for (dis in 1:nbD){
      AIC[dis]<- -2*(llk.list[[dis]]-numb.par)
      BIC[dis]<- -2*llk.list[[dis]]+numb.par*log(nbT)}
    
    ## 3.Stop iteration          ----------------------------
    New.param <-  c(esLambda, transPr)
    crit <- New.param - Old.param
    Old.param  <- New.param
    
    if (iter > 1 && max(abs(crit)) <= threshold) break
  }
  return(list(lpostPr = lpostPr.list, initPr = initPr,  transPr = transPr,  esLambda = esLambda,
              iterstop = iter,crit=crit, AIC=AIC, BIC=BIC))
}


mod<- CHMM_VEM(X=X, nb.states=3, itmax = 100, threshold = 0.00001, nbD=nbD, omega=omega)
mod