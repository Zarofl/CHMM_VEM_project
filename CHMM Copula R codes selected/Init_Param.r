#necessary libraries
library(Rmpfr)
library(dplyr)

### --Function to generate initial parameters of CopulaCHMM

init.VEM <- function(X, nb.states, nbD, nbT){
  #emission prob-------------- 
  #Will hold emission parameters per feature and state
  #rows are diseases, columns are states
  esLambda<-matrix(NA, nbD, nb.states)
  #make sure all the intial emission start from the same mean
  #Across all states
  for (i in 1:nbD){
    esLambda[i,]<-mean(X[,i])}
  
  #Creating the estimated next state of emission = prev+50%
  for (j in 1:(nb.states-1)){
    esLambda[,j+1]<-esLambda[,j]+0.5*esLambda[,j]}

  #set.seed(123)

  # Transition init: high self-probability
  mat.tmp <- matrix(runif(nb.states^2), ncol = nb.states) + diag(rep(50, nb.states))
  transPr <- mat.tmp / rowSums(mat.tmp)

  # initial distribution: stationary of transPr, then clamp-----------------------------------------
  eigenvalues <- round(eigen(t(transPr))$values, 3)
  pos  <- which(eigenvalues == 1.000)
  nuHMM <- eigen(t(transPr))$vectors[, pos]
  initPr <- pmax(as.numeric(nuHMM / sum(nuHMM)), 0)
  initPr <- initPr / sum(initPr)
  
  # postPr (random simplex at each t)  ------------------------------------------------------
  postPr <- list()
  for(dis in 1:nbD){
    #tau.tmp <- data.frame(matrix(runif(nbT * nb.states, 0, 1), nbT, nb.states))
    #set.seed(144)
    tau.tmp <- matrix(runif(nbT * nb.states, 0, 1), nbT, nb.states)
    tau.tmp <- tau.tmp / rowSums(tau.tmp)
    postPr[[dis]] <- tau.tmp
  }
  return(list(esLambda = esLambda, transPr = transPr, initPr = initPr,
              postPr = postPr))
  
}
