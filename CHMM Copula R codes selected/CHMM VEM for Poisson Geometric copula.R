#necessary libraries 
library(Rmpfr)
library(dplyr)

##--Function to generate initial parameters of CopulaCHMM
# init.VEM <- function(X, nb.states, nbD, nbT){
#   #emission prob--------------
#   esLambda<-matrix(NA,nbD,nb.states)  #rows are diseases, columns are states
#   for (i in 1:nbD){
#     esLambda[i,]<-mean(X[,i])
#   }
  
#   for (j in 1:(nb.states-1)){
#     esLambda[,j+1]<-esLambda[,j]+0.5*esLambda[,j]
#   }
  
  
#   #set.seed(123)
#   mat.tmp <- matrix(runif(nb.states^2), ncol = nb.states) + diag(rep(50, nb.states))
#   transPr <- mat.tmp / rowSums(mat.tmp)
#   # initial distribution -----------------------------------------
#   eigenvalues <- round(eigen(t(transPr))$values, 3)
#   pos  <- which(eigenvalues == 1.000)
#   nuHMM <- eigen(t(transPr))$vectors[, pos] 
#   initPr <- pmax(as.numeric(nuHMM / sum(nuHMM)), 0)
#   initPr <- initPr / sum(initPr)
  
#   # postPr  ------------------------------------------------------
#   postPr <- list()
#   for(dis in 1:nbD){
#     #tau.tmp <- data.frame(matrix(runif(nbT * nb.states, 0, 1), nbT, nb.states))
#     #set.seed(144)
#     tau.tmp <- matrix(runif(nbT * nb.states, 0, 1), nbT, nb.states)
#     tau.tmp <- tau.tmp / rowSums(tau.tmp)
#     postPr[[dis]] <- tau.tmp
#   }
#   return(list(esLambda = esLambda, transPr = transPr, initPr = initPr,
#               postPr = postPr))
  
# }

#Function to calculate log-forward probabilities
ForwardR <- function(initPr,emisPrW, transPr, nbT, nb.states) {
  lFpr<-matrix(NA,nbT,nb.states)
  
  for (i in 1:nb.states){
    print(initPr[i])
    if (initPr[i]==0) {
      initPr[i]=initPr[i]+0.0001
    } 
    if (initPr[i]==1){
      initPr[i]=initPr[i]-0.0002}
  }
  
  foo<-initPr
  sumfoo<-sum(foo)
  lscale<-log(sumfoo)
  foo<-foo/sumfoo
  lFpr[1,]<-log(foo)+lscale
  for (t in 2:nbT){
    z1<-colSums(foo%*%transPr)
    foo<-emisPrW[t,]*z1
    sumfoo<-sum(foo)
    lscale<-lscale+log(sumfoo)
    foo<-foo/sumfoo
    lFpr[t,]<-log(foo)+lscale
  }
  return(list(lFpr = lFpr))
}


##Function to calculate log-backward probabilities

BackwardR <- function(lFpr, transPr, nbT, nb.states) {
  
  #lFpr<-resF$lFpr
  
  lpostPr<-matrix(NA,nbT,nb.states)
  
  lGpr<-matrix(NA,nbT,nb.states)
  lGp<-matrix(NA,nbT,nb.states)
  lpostPr[nbT,]<-lFpr[nbT,]
  
  
  
  for (t in (nbT-1):1){
    #lGp[t,]<-lFpr[t,]+log(transPr)
    c<-max(lFpr[t,]+log(transPr))
    #print(c)
    
    d<-asNumeric(exp(mpfr(lFpr[t,]+log(transPr)-c,5)))  ##here we use mpfr package function to precisely calculate exp function
    #d<-matrix(d,nbT,nb.states)
    
    lGpr[t+1,]<-c+log(colSums(d)) #think
    #print(lGpr[t+1,])
    
    #Gpr<-exp(lGpr[t+1,])
    #print(exp(lFpr[t,]))
    #print(lFpr[t,])
    #delta  <- exp(lFpr[t,])*transPr/Gpr*lpostPr[t+1,]
    
    ldelta  <- lFpr[t,]+log(transPr)-lGpr[t+1,]+lpostPr[t+1,] 
    #print(ldelta)
    c2<-max(ldelta)
    
    d2<-asNumeric(exp(mpfr(ldelta-c2,5)))
    #d2<-matrix(d2,nb.states,nb.states)
    
    
    
    lpostPr[t,]<-c2+log(rowSums(d2))
    #print(ldelta)
    #delta<-exp(ldelta)
    #print(delta)
    #lpostPr[t,]<-rowSums(delta)
    #lpostPr[t,]<-log(lpostPr[t,])
    
    
  }
  return(list(lGpr=lGpr,lpostPr=lpostPr))
  
}
########################################################
##Our main function of Copula CHMM, includes EM algorithm, calculates AIC and BIC, 
#also calculates RSS (however I will change it for Poisson, so for now we have NA )

CHMM_VEM <- function(X, nb.states, 
                     itmax = 5e2,
                     threshold = 1e-7,nbD,omega, res.init){ 
  nbT <- nrow(X) 


# pick eigenvector corresponding to eigenvalue closest to 1
  # eig <- eigen(t(res.init$transPr))
  # pos <- which.min(abs(eig$values - 1))
  # nuHMM <- Re(eig$vectors[, pos])
  # res.init$initPr <- pmax(as.numeric(nuHMM), 0)
  # res.init$initPr <- res.init$initPr / sum(res.init$initPr)

  
  #initial parameters
  #res.init <- init.VEM(X=X, nb.states=nb.states, nbD=nbD, nbT=nbT)

  print(res.init)
  esLambda <- res.init$esLambda
  print(esLambda)
  transPr <- res.init$transPr
  print(transPr)
  initPr <- res.init$initPr
  print(initPr)
  #postPr.last <- res.init$postPr
  #print(postPr.last)
  
  
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
    
    
    #-----Calculate geometric copula
    copula<-matrix(0,nbD,nbD)
    prodcopula<-rep(NA,nbD)
    RSS<-matrix(NA,nbD,nbD)
    for (dis in 1:nbD) {
      for (dis2 in c(1:nbD)[-dis]){ 
        om<-omega[dis,dis2]
        geometric_copula=matrix(NA,nb.states,nb.states) #n is the number of trials
        geometric_copula[1,1]<- 1/3*2*om/(2*om+sqrt(8*om+1)+1)
        geometric_copula[1,2]<- 1/3*(sqrt(8*om+1)+1)/2/(2*om+sqrt(8*om+1)+1)
        geometric_copula[1,3]<- 1/3*(sqrt(8*om+1)+1)/2/(2*om+sqrt(8*om+1)+1)
        geometric_copula[2,1]<- 1/3*(sqrt(8*om+1)+1)/2/(2*om+sqrt(8*om+1)+1)
        geometric_copula[2,2]<- 1/3*sqrt(om)*(sqrt(8*om+1)+1)^2/4/(sqrt(om)+1)/(2*om+sqrt(8*om+1)+1)
        geometric_copula[2,3]<- 1/3*(4*om-1-sqrt(8*om+1))/4/(sqrt(om)-1)/(sqrt(om)+1)^2
        geometric_copula[3,1]<- 1/3*(sqrt(8*om+1)+1)/2/(2*om+sqrt(8*om+1)+1)
        geometric_copula[3,2]<- 1/3*(4*om-1-sqrt(8*om+1))/4/(sqrt(om)-1)/(sqrt(om)+1)^2
        geometric_copula[3,3]<- 1/3*sqrt(om)*(sqrt(8*om+1)+1)^2/4/(sqrt(om)+1)/(2*om+sqrt(8*om+1)+1)
        
        
        #copula[dis,dis2]<- geometric_copula[1,1]
        #copula[dis,dis2]<- geometric_copula[1,2]
        #copula[dis,dis2]<- geometric_copula[2,2]
        #copula[dis,dis2]<- geometric_copula[1,3]
        copula[dis,dis2]<- max(geometric_copula)
        
      }
      
      #prodcopula[dis]<-prod(copula[dis,-dis]) 
      prodcopula[dis]<-max(geometric_copula)
      
      lpost_z<-rep(NA,nbT)
      
      trsTmp <- matrix(0, nb.states, nb.states) 
      
      #Calculate emission probabilitie that follow Poisson

      print(esLambda)
      #
      for (t in 1:nbT)  {
        emisPr[t,] <- dpois(X[t,dis],lambda=esLambda[dis,]) #for each disease 
      } 
      
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
        a1<-exp(mpfr(lpostPr.tmp[t,],5))
        a1<-asNumeric(a1)
        
        lpost_z[t]<-max(lpostPr.tmp[t,])+log(sum(a1-max(lpostPr.tmp[t,]))) 
      }  
      
      
      lpostPr <- lpostPr.tmp - lpost_z 
      
      
      # postPr <- postPr.tmp / rowSums(postPr.tmp)  
      
      #lpostPr <- lpostPr.tmp / rowSums(lpostPr.tmp)
      
      lpostPr.list[[dis]] <- lpostPr
      lFpr <- resF$lFpr
      lGpr <- resB$lGpr
      
      #--Here we calculate matrix using formulas in backward recursion, actually we calculate conditional expectation
      # c<-asNumeric(exp(mpfr(lFpr[-nbT,],5)))
      c<-exp(mpfr(lFpr[-nbT,],5))
      
      #  d<-asNumeric(exp(mpfr(lpostPr[-1,] - lGpr[-1,],5)))
      d<-exp(mpfr(lpostPr[-1,] - lGpr[-1,],5))  
      #d<-asNumeric(d)
      
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
    
    
    #  order the parameters          ----------------------------
    #ordLambda <- order(colSums(esLambda))
    #ordLambda <- order(esLambda[1,])
    #esLambda<- esLambda[,ordLambda]
    #transPr <- transPr[ordLambda,ordLambda]
    
    ## Calculate criteria    RSS     -------------------------- I will update it 
    
    for (dis in 1:nbD) {
      for (dis2 in c(1:nbD)[-dis]){
        a<-as.numeric(exp(mpfr(lpostPr.list[[dis]],5)))
        a<-matrix(a,nbT,nb.states)
        b<-as.numeric(exp(mpfr(lpostPr.list[[dis2]],5)))
        b<-matrix(b,nbT,nb.states)
        RSS[dis,dis2]<-sum(a*(X[,dis]-mean(esLambda[dis,]))^2)+sum(b*(X[,dis2]-mean(esLambda[dis2,]))^2)
        #RSS[dis,dis2]<-sum(exp(mpfr(lpostPr.list[[dis]],5))*(X[,dis]-nbI[dis]*esProb[dis,])^2)+sum(exp(mpfr(lpostPr.list[[dis2]],5))*(X[,dis2]-nbI[dis2]*esProb[dis2,])^2)
      }
    }
    
    ## Calculate AIC and BIC
    
    for (dis in 1:nbD){
      c<-max(resF.list[[dis]][nbT,])
      llk.list[[dis]]<-c+log(sum(exp(resF.list[[dis]][nbT,]-c)))
    }
    
    numb.par<-nbD*(nbD-1)/2+nbD*nb.states+nb.states*(nb.states-1)+nb.states-1 
    
    for (dis in 1:nbD){
      AIC[dis]<- -2*(llk.list[[dis]]-numb.par)
      BIC[dis]<- -2*llk.list[[dis]]+numb.par*log(nbT)
    }
    
    
    ## 3.Stop iteration          ----------------------------
    New.param <-  c(esLambda, transPr)
    crit <- New.param - Old.param
    Old.param  <- New.param
    
    if (iter > 1 && max(abs(crit)) <= threshold) break()
  }
  return(list(lpostPr = lpostPr.list, initPr = initPr,  transPr = transPr,  esLambda = esLambda,
              iterstop = iter,crit=crit, AIC=AIC, BIC=BIC))
}

# mod<- CHMM_VEM(X=X, nb.states=3, itmax = 100, threshold = 0.00001, nbD=nbD, omega=omega)
# mod

#######################
print(summary(X))

#omega<-matrix(c(NA,1.5,1.5,NA),nrow=nbD,ncol=nbD)

res.init<-list()
res.init<-list(matrix(c(14,23,30,14,24,30),2,3,byrow=T),matrix(c(0.8,0.1,0.1,0.1,0.8,0.1,0.1,0.1,0.8),byrow=T,3,3),
               c(0.8,0.1,0.1))
names(res.init) <- c("esLambda", "transPr", "initPr")

mod<- CHMM_VEM(X=X, nb.states=3, itmax = 100, threshold = 0.00001, nbD=nbD, omega=omega, res.init=res.init)
mod
