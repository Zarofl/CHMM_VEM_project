#necessary libraries
library(Rmpfr)
library(dplyr)

#Function to calculate log-forward probabilities
ForwardR <- function(initPr,emisPrW, transPr, nbT, nb.states) {
  lFpr<-matrix(NA,nbT,nb.states)
  
  for (i in 1:nb.states){
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