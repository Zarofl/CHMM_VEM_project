#############################################################################################
##########------Hocam this code reads real data and ready to use for model-------------------------------

###---Data
Y<-read.csv("D://Copula - QU R code/CHMM_VEM_project/CHMM Copula R codes selected/data.csv",header=T,sep=",")
str(Y)
summary(Y)

NR<-nrow(Y)
NC<-ncol(Y)

#This matrix is suitable for the dataset to keep outliers 
#It's adviced to keep the values between 0-100. 
Y2<-matrix(NA,NR,NC)
for (j in 2:NC){
  for (i in 2:NR){
    Y2[i,j]<-Y[i,j]-Y[i-1,j]
  }
}
Y2<-as.data.frame(Y2)
Y2[1,]<-Y[1,]
Y2[,1]<-Y[,1]
Y2<-Y2[,-5]

hist(Y2[,2])
hist(Y2[,3])
hist(Y2[,4])

Y2[,2]<-round(Y2[,2]/100,0)
Y2[,3]<-round(Y2[,3]/100,0)
Y2[,4]<-round(Y2[,4]/100,0)
#Y2[,4]<-rpois(72,49)


nrow(Y2)
nbD<-ncol(Y2)-1

summary(Y2)

X<-Y2
X<-X[,-1]
#X<-X[,-3]   #without 3rd disease


hist(X$V2)
hist(X$V3)
hist(X$V4)

###Start the function
nrow(X)
nbD<-ncol(X)

# nbI<-rep(NA,nbD)
# for (i in 1:nbD){
#   nbI[i]<-max(X[,i])+0.01*sd(X[,i])
# }
# nbI<-round(nbI,0)

omega<-matrix(c(NA,0.8,1.2,0.8,NA,0.7,1.2,0.7,NA),nrow=nbD,ncol=nbD)
#omega<-matrix(c(NA,0.8,0.8,NA),nrow=nbD,ncol=nbD)

