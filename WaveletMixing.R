# Oliver Orejola
# Created: 7/2/21
# Updated: 9/16/21
# Wavelet ESD for Mixed Hurst Parameters



#########################[ MIXED ]#############################
#######################[ Run First ]###########################
library(somebm)
library(wavelets)
library(diptest)
library(ggplot2)

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



N_mixed_wavelet_sampe_covariance <- function(level,scale,p,Hursts){
  N = 2^level
  n_j = N/(2^scale) #effective sample size
  X_MIXED <- c()
  for(i in 1:p){
    MIXED<- c()
    H <- sample(Hursts,1)
    MIXED<-  N^H*fbm(hurst = H, n = N)
    wave_MIXED <- dwt(MIXED, filter = "d4", n.levels = level, boundary = "reflection")
    X_MIXED<- append(X_MIXED,unlist(wave_MIXED@W[scale]))
  }
  Xmat_MIXED <- matrix( X_MIXED, n_j, p,byrow = FALSE)  
  W_MIXED <-  crossprod(Xmat_MIXED)
  return(log(eigen(1/n_j*W_MIXED)$values))
}


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

monte_carlo_power_N_hursts <- function(runs,significance,level,scale,p,hursts){
  count = 0
  for(i in 1:runs){
    E_MIXED <- N_mixed_wavelet_sampe_covariance(level,scale,p,hursts)
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


Etest <- N_mixed_wavelet_sampe_covariance(16,8,50,c(0.2,0.8))
hist(Etest)


m <- min(Etest)
M <- max(Etest)

hist(Etest,probability=TRUE
     ,breaks = seq(m-1,M+1, by = 0.25)
     ,xlim = c(m,M+5)
)

monte_carlo_power_N_hursts(10,0.05,15,7,50,c(0.2,0.8,0.4,0.6))

E_MIXED <- mixed_wavelet_sampe_covariance(level,scale,p,H0,H1,mixingrate)
m <- min(E_MIXED)
M <- max(E_MIXED)

hist(E_MIXED,probability=TRUE
     ,breaks = seq(m-1,M+1, by = 0.25)
     ,xlim = c(m,M+5)
     ,xlab = "Eigenvalues"
     ,main = paste("H_0=",H0,", H_1=",H1, ", n_j=",n_j,", p/n_j=",alpha)
)

#-------------------------------------- dimension = 32
#----------------------- 2 hurst parameters 
for(i in 1:6){
  level = 10+i
  scale = 4+i
  N = 2^level #Path size
  n_j = N/(2^scale)#effective sample size
  alpha = 0.5 # effective sample size and dimension ratio, must be greater than 1 for rank sufficiency 
  p = floor(alpha*n_j) #
  power <-monte_carlo_power_N_hursts(1000,0.05,level,scale,p,c(0.27,0.75))
  print(paste(scale,level,p,power))
}
#----------------------- 4 hurst parameters 
for(i in 1:6){
  level = 10+i
  scale = 4+i
  N = 2^level #Path size
  n_j = N/(2^scale)#effective sample size
  alpha = 0.5 # effective sample size and dimension ratio, must be greater than 1 for rank sufficiency 
  p = floor(alpha*n_j) #
  power <-monte_carlo_power_N_hursts(1000,0.05,level,scale,p,c(0.2,0.4,0.6,0.8))
  print(paste(scale,level,p,power))
}

#----------------- 1 hurst paramter
for(level in 10:16){
  scale = level-6
  N = 2^level #Path size
  n_j = N/(2^scale)#effective sample size
  alpha = 0.5 # effective sample size and dimension ratio, must be greater than 1 for rank sufficiency 
  p = floor(alpha*n_j) #
  power <-monte_carlo_power_N_hursts(1000,0.05,level,scale,p,c(0.5))
  print(paste(scale,level,p,power))
}


#----------------- 3 hurst paramter
for(level in 10:16){
  scale = level-6
  N = 2^level #Path size
  n_j = N/(2^scale)#effective sample size
  alpha = 0.5 # effective sample size and dimension ratio, must be greater than 1 for rank sufficiency 
  p = floor(alpha*n_j) #
  power <-monte_carlo_power_N_hursts(1000,0.05,level,scale,p,c(0.2,0.5,0.8))
  print(paste(scale,level,p,power))
}

#----------------- 6 hurst paramter
for(level in 10:16){
  scale = level-6
  N = 2^level #Path size
  n_j = N/(2^scale)#effective sample size
  alpha = 0.5 # effective sample size and dimension ratio, must be greater than 1 for rank sufficiency 
  p = floor(alpha*n_j) #
  power <-monte_carlo_power_N_hursts(10,0.05,level,scale,p,c(0.1,0.25,0.40,0.55,0.70,0.85))
  print(paste(scale,level,p,power))
}


par(mfrow = c(1,1))
E_MIXED <- N_mixed_wavelet_sampe_covariance(15,6,128,c(0.1,0.25,0.40,0.55,0.70,0.85))
m <- min(E_MIXED)
M <- max(E_MIXED)
N = 2^15 #Path size
n_j = N/(2^6)
pvalue = dip.test(E_MIXED)$p.value
hist(E_MIXED,probability=TRUE
     ,breaks = seq(m-1,M+1, by = 0.20)
     ,xlim = c(m,M)
     ,xlab = "Eigenvalues"
     ,main = paste("p-value",pvalue)
)


#--------- Plots -----------------------
par(mfrow = c(2,2))
xaxis<- c(10:16)
#--------- 1 hurst ---------------------
power_1_hurst <- c(0,0,0,0,0,0,0)
plot(xaxis,power_1_hurst,type="o"
     ,xlab = "level"
     ,ylab = "Proportion of Rejections",ylim=c(0,1)
     ,main = "One Hurst Parameter")

#--------- 2 hurst ---------------------
power_2_hurst <- c(0,0.068,0.603,0.844,0.931,0.943,0.973)
plot(xaxis,power_2_hurst,type="o"
     ,xlab = "level"
     ,ylab = "Proportion of Rejections",ylim=c(0,1)
     ,main = "Two Hurst Parameters")

#--------- 3 hurst ---------------------
power_3_hurst <- c(0,0.007,0.039,0.142,0.312,0.551,0.734)
plot(xaxis,power_3_hurst,type="o"
     ,xlab = "level"
     ,ylab = "Proportion of Rejections",ylim=c(0,1)
     ,main = "Three Hurst Parameters")

#--------- 4 hurst ---------------------
power_4_hurst <- c(0,0,0.006,0.023,0.052,0.104,0.146)
plot(xaxis,power_4_hurst,type="o"
     ,xlab = "level"
     ,ylab = "Proportion of Rejections",ylim=c(0,1)
     ,main = "Four Hurst Parameters")





#--------- 6 hurst ---------------------

power_6_hurst <- c(0,0,0,0,0.052,0.104,0.2)
plot(xaxis,power_6_hurst,type="o"
     ,xlab = "level",ylim=c(0,1)
     ,ylab = "Proportion of Rejections"
     ,main = "6 Hurst Parameters")

