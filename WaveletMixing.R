# Oliver Orejola
# Created: 7/2/21
# Updated: 7/18/21
# Wavelet ESD for Mixed Hurst Parameters



#########################[ MIXED ]#############################
#######################[ Run First ]###########################
library(somebm)
library(wavelets)


H0 = 0.1
H1 = 0.25
level = 10
scale = 6
mixingrate = 0.5

path_size = 2^level
N = path_size 
n_j = N/(2^scale)
X_MIXED <- c()
X_0 <- c()
X_1 <- c()
for(i in 1:n_j){
  fbmSIM_0 <- N^H0*fbm(hurst = H0, n = N)
  fbmSIM_1 <- N^H1*fbm(hurst = H1, n = N)
  MIXED<- c()
  if( rbinom(1,1,mixingrate) == 1){
    MIXED<- append(MIXED, fbmSIM_1)
  }
  else{
    MIXED<- append(MIXED, fbmSIM_0)
  }
  wave_MIXED <- dwt(MIXED, filter = "d4", n.levels = level, boundary = "reflection")
  wave_0 <- dwt(fbmSIM_0, filter = "d4", n.levels = level, boundary = "reflection")
  wave_1 <- dwt(fbmSIM_1, filter = "d4", n.levels = level, boundary = "reflection")
  X_MIXED<- append(X_MIXED,unlist(wave_MIXED@W[scale]))
  X_0 <- append(X_0,unlist(wave_0@W[scale]))
  X_1<- append(X_1,unlist(wave_1@W[scale]))
}


Xmat_MIXED <- matrix( X_MIXED, n_j, n_j,byrow = TRUE)  
W_MIXED <-  crossprod(Xmat_MIXED)
E_MIXED <-log(eigen(1/n_j*W_MIXED)$values)

Xmat_0 <- matrix( X_0, n_j, n_j,byrow = TRUE)  
W_0 <-  crossprod(Xmat_0)
E_0 <-log(eigen(1/n_j*W_0)$values)

Xmat_1 <- matrix( X_1, n_j, n_j,byrow = TRUE)  
W_1 <-  crossprod(Xmat_1)
E_1 <-log(eigen(1/n_j*W_1)$values)


#########################
par(mfrow = c(4,3))

m <- min(E_MIXED,E_0,E_1)
M <- max(E_MIXED,E_0,E_1)

hist(E_MIXED,probability=TRUE
     ,breaks = seq(m-1,M+1, by = 0.75)
     ,xlim = c(m,M+5)
     ,main = paste("MIXED, n=2^",level,", j=",scale)
)
hist(E_0 ,probability=TRUE
     ,breaks = seq(m,M+1, by = 0.75)
     ,xlim = c(m,M+5)
     ,main = paste("H_0=",H0,", n=2^",level,", j=",scale)
)
hist(E_1 ,probability=TRUE
     ,breaks = seq(m-1,M+1, by = 0.75)
     ,xlim = c(m,M+5)
     ,main = paste("H_1=",H1,", n=2^",level,", j=",scale)
)
