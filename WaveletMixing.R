# Oliver Orejola
# Created: 7/2/21
# Updated: 9/16/21
# Wavelet ESD for Mixed Hurst Parameters



#########################[ MIXED ]#############################
#######################[ Run First ]###########################
library(somebm)
library(wavelets)
library(diptest)

H0 = 0.25
H1 = 0.75
level = 14
scale = 8
mixingrate = 0.5

N = 2^level #Path size
n_j = N/(2^scale)#effective sample size

alpha = 0.5 # effective sample size and dimension ratio, must be greater than 1 for rank sufficiency 
p = floor(alpha*n_j) #dimension i.e. number of fBM we sample


mixed_wavelet_sampe_covariance <- function(level,scale,p,H0,H1,mixingrate){
  N = 2^level
  n_j = N/(2^scale) #effective sample size
  X_MIXED <- c()
  for(i in 1:p){
    MIXED<- c()
    if( rbinom(1,1,mixingrate) == 1){
      fbmSIM_1 <- N^H1*fbm(hurst = H1, n = N)
      MIXED<- fbmSIM_1
    }
    else{
      fbmSIM_0 <- N^H0*fbm(hurst = H0, n = N)
      MIXED<- fbmSIM_0
    }
    wave_MIXED <- dwt(MIXED, filter = "d4", n.levels = level, boundary = "reflection")
    X_MIXED<- append(X_MIXED,unlist(wave_MIXED@W[scale]))
  }
  Xmat_MIXED <- matrix( X_MIXED, n_j, p,byrow = FALSE)  
  W_MIXED <-  crossprod(Xmat_MIXED)
  return(log(eigen(1/n_j*W_MIXED)$values))
}

E_MIXED <- mixed_wavelet_sampe_covariance(level,scale,p,H0,H1,mixingrate)

monte_carlo_power_two_hurst <- function(runs,significance,level,scale,p,mixingrate){
  count = 0
  for(i in 1:runs){
    E_MIXED <- mixed_wavelet_sampe_covariance(level,scale,p,H0,H1,mixingrate)
    pvalue<-dip.test(E_MIXED)$p.value
    if( pvalue<= significance){
      count = count + 1
    }
  }
  return(count/runs)
}


#-----------------------------------------------------

H0 = 0.25
H1 = 0.75
mixingrate = 0.5

for(i in 1:5){
  level = 9+i
  scale = 1+i
  N = 2^level #Path size
  n_j = N/(2^scale)#effective sample size
  alpha = 0.5 # effective sample size and dimension ratio, must be greater than 1 for rank sufficiency 
  p = floor(alpha*n_j) #
  power <-monte_carlo_power_two_hurst(1000,0.05,level,scale,p,0.5)
  print(paste(scale,level,p,power))
}

for(i in 1:5){
  level = 9+i
  scale = 2+i
  N = 2^level #Path size
  n_j = N/(2^scale)#effective sample size
  alpha = 0.5 # effective sample size and dimension ratio, must be greater than 1 for rank sufficiency 
  p = floor(alpha*n_j) #
  power <-monte_carlo_power_two_hurst(1000,0.05,level,scale,p,0.5)
  print(paste(scale,level,p,power))
}


plot(c(10,11,12,13,14),c(0,0,0.776,0.988,0.997),type = "o")




pvalue<-dip.test(E_MIXED)$p.value

alpha = 0.5
significance = 0.05 #Rejection  region

for(i in 9:12){
  for(j in 4:(i-1)){
    p = floor( alpha * 2^(i-j))
    E_MIXED <- mixed_wavelet_sampe_covariance(i,j,p,H0,H1,mixingrate)
    pvalue<-dip.test(E_MIXED)$p.value
    
    if( pvalue<= significance){
      test = "Reject"
    }
    else{
      test = "Accept"
    }
    print(paste(round(j,digits=1),as.integer(i),as.integer(p), pvalue, test))
  }
}
#########################
par(mfrow = c(1,1))



m <- min(E_MIXED)
M <- max(E_MIXED)

hist(E_MIXED,probability=TRUE
     ,breaks = seq(m-1,M+1, by = 0.25)
     ,xlim = c(m,M+5)
     ,xlab = "Eigenvalues"
     ,main = paste("H_0=",H0,", H_1=",H1, ", n_j=",n_j,", p/n_j=",alpha)
)




