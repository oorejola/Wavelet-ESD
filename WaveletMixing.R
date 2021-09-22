# Oliver Orejola
# Created: 7/2/21
# Updated: 9/21/21
# Wavelet ESD for Mixed Hurst Parameters
#------------------------- [ Libraries ]----------------------------
library(somebm)
library(wavelets)
library(diptest)
library(ggplot2)


#------------------------- [ Functions ]----------------------------

coordinate_matrix <- function(n){
  mat <- c()
  for(i in 1:n){
    x <- rnorm(n)
    col <- x/sqrt(sum(x^2)) 
    mat <- append(mat,col)
  }
  return(matrix(mat, n, n,byrow = FALSE))
}

wavelet_sampe_covariance<- function(level,scale,p,Hursts, mixed = FALSE){
  N = 2^level
  n_j = N/(2^scale) #effective sample size
  wavelet_list <- c() 
  for(i in 1:p){
    H <- sample(Hursts,1)
    #MIXED<-  N^H*fbm(hurst = H, n = N) 
    fractional_brownian <-  N^H*fbm(hurst = H, n = N) 
    wavelet_transform <- dwt(fractional_brownian, filter = "d4", n.levels = level, boundary = "reflection")
    wavelet_list <- append(wavelet_list,unlist(wavelet_transform@W[scale]))
  }
  wavelet_mat<- matrix(wavelet_list, n_j, p,byrow = FALSE)
  
  if(mixed == TRUE){ 
    P <- coordinate_matrix(p)
    mixed_wavlet_mat <- wavelet_mat %*% t(P)
    mixed_wavlet_product <- crossprod(mixed_wavlet_mat)
    eigenvalues <- eigen(1/n_j*mixed_wavlet_product)$values
  }
  else{
    wavelet_product <-  crossprod(wavelet_mat)
    eigenvalues <- eigen(1/n_j*wavelet_product)$values
  }
  return(eigenvalues)
}

monte_carlo_probability_of_rejection <- function(runs,significance,level,scale,p,hursts){
  count = 0
  for(i in 1:runs){
    E_MIXED <- log(wavelet_sampe_covariance(level,scale,p,hursts))
    pvalue<-dip.test(E_MIXED)$p.value
    if( pvalue<= significance){
      count = count + 1
    }
  }
  proportion_of_rejections = count/runs
  return(proportion_of_rejections)
}


#-----------------------------------------------------
H0 = 0.25
H1 = 0.75
level = 14
scale = 8
mixingrate = 0.5

N = 2^level #Path size
n_j = N/(2^scale)#effective sample size

alpha = 0.5 # effective sample size and dimension ratio, must be greater than 1 for rank sufficiency 
p = floor(alpha*n_j) #dimension i.e. number of fBM we sample

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
  power <-monte_carlo_probability_of_rejection(1000,0.9,level,scale,p,c(0.25,0.75))
  print(paste(scale,level,p,power))
}



plot(c(10,11,12,13,14),c(0,0,0.776,0.988,0.997),type = "o")


pvalue<-dip.test(E_MIXED)$p.value

alpha = 0.5
significance = 0.05 #Rejection  region

for(i in 9:16){
  for(j in 4:(i-1)){
    p = floor( alpha * 2^(i-j))
    #E_MIXED <- log(wavelet_sampe_covariance(i,j,p,H0,H1,mixingrate))
    #pvalue<-dip.test(E_MIXED)$p.value
    pvalue = 0
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


Etest <- log(wavelet_sampe_covariance(16,8,50,c(0.2,0.8)))
hist(Etest)


m <- min(Etest)
M <- max(Etest)

hist(Etest,probability=TRUE
     ,breaks = seq(m-1,M+1, by = 0.25)
     ,xlim = c(m,M)
)

monte_carlo_probability_of_rejection(10,0.05,15,7,50,c(0.2,0.8,0.4,0.6))



#-------------------------------------- dimension = 32
#----------------------- 2 hurst parameters 
for(level in 10:13){
  scale = level - 6
  N = 2^level #Path size
  n_j = N/(2^scale)#effective sample size
  alpha = 0.5 # effective sample size and dimension ratio, must be greater than 1 for rank sufficiency 
  p = floor(alpha*n_j) #
  power <-monte_carlo_probability_of_rejection(100,0.90,level,scale,p,c(0.27,0.75))
  print(paste(scale,level,p,power))
}
#---- FIXED SCALE

#---- Triple Limit NULL
for(scale in 5:8){
  level = scale*2
  p = 2^(level- 1-scale) #given alpha = 1/2
  power <-monte_carlo_probability_of_rejection(1000,0.90,level,scale,p,c(0.5))
  print(paste(scale,level,p,power))
}

#---- Triple Limit 2 Hurst
for(scale in 5:8){
  level = scale*2
  p = 2^(level- 1-scale) #given alpha = 1/2
  power <-monte_carlo_probability_of_rejection(10,0.90,level,scale,p,c(0.27,0.75))
  print(paste(scale,level,p,power))
}
#---- Triple Limit 4 Hurst
for(scale in 5:8){
  level = scale*2
  p = 2^(level- 1-scale) #given alpha = 1/2
  power <-monte_carlo_probability_of_rejection(10,0.90,level,scale,p,c(0.2,0.4,0.6,0.8))
  print(paste(scale,level,p,power))
}


#----------------------- 4 hurst parameters 
for(level in 10:13){
  scale = level - 6
  N = 2^level #Path size
  n_j = N/(2^scale)#effective sample size
  alpha = 0.5 # effective sample size and dimension ratio, must be greater than 1 for rank sufficiency 
  p = floor(alpha*n_j) #
  power <-monte_carlo_probability_of_rejection(100,0.90,level,scale,p,c(0.2,0.4,0.6,0.8))
  print(paste(scale,level,p,power))
}

#----------------- 1 hurst paramter (NULL)
for(level in 10:13){
  scale = level-6
  N = 2^level #Path size
  n_j = N/(2^scale)#effective sample size
  alpha = 0.5 # effective sample size and dimension ratio, must be greater than 1 for rank sufficiency 
  p = floor(alpha*n_j) #
  power <-monte_carlo_probability_of_rejection(100,0.95,level,scale,p,c(0.5))
  print(paste(scale,level,p,power))
}


#----------------- 3 hurst paramter
for(level in 10:16){
  scale = level-6
  N = 2^level #Path size
  n_j = N/(2^scale)#effective sample size
  alpha = 0.5 # effective sample size and dimension ratio, must be greater than 1 for rank sufficiency 
  p = floor(alpha*n_j) #
  power <-monte_carlo_probability_of_rejection(1000,0.05,level,scale,p,c(0.2,0.5,0.8))
  print(paste(scale,level,p,power))
}

#----------------- 6 hurst paramter
for(level in 10:16){
  scale = level-6
  N = 2^level #Path size
  n_j = N/(2^scale)#effective sample size
  alpha = 0.5 # effective sample size and dimension ratio, must be greater than 1 for rank sufficiency 
  p = floor(alpha*n_j) #
  power <-monte_carlo_probability_of_rejection(100,0.05,level,scale,p,c(0.1,0.25,0.40,0.55,0.70,0.85))
  print(paste(scale,level,p,power))
}


par(mfrow = c(1,1))
E_MIXED <- log(wavelet_sampe_covariance(15,6,128,c(0.1,0.25,0.40,0.55,0.70,0.85)))
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


#---------[ Plots ] (Fixed dimension)-----------------------
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




E<- log(wavelet_sampe_covariance(16,8,64,c(0.5)))/(j*log(2))
m=min(E)
M=max(E)
hist(E,probability=TRUE
     ,breaks = seq(m-1,M+1, by = 0.25)
#     ,xlim = c(m,M)
     ,xlab = "Eigenvalues"
)
mean(E)

#-------- Investigate Scales Triple Limit---------------------
for(scale in 5:8){
  level = scale*2
  p = 2^(level- 1-scale)
  E <- log(wavelet_sampe_covariance(level,scale,p,c(0.5)))/(scale*log(2))
  mean <- mean(E)
  print(paste(scale,level,as.integer(p), mean))
}


for(scale in 5:8){
  level = scale*2
  p = 2^(level- 1-scale)
  E <- log(wavelet_sampe_covariance(level,scale,p,c(0.75)))/(scale*log(2))
  mean <- mean(E)
  print(paste(scale,level,as.integer(p), mean))
}

for(scale in 5:13){
  level = 18
  p = 32
  E <- log(wavelet_sampe_covariance(level,scale,p,c(0.75)))/(scale*log(2))
  mean <- mean(E)
  print(paste(scale,level,as.integer(p), mean))
}
for(scale in 5:13){
  level = 18
  p = 32
  E <- log(wavelet_sampe_covariance(level,scale,p,c(0.75)))/(scale*log(2))
  mean <- mean(E)
  print(paste(scale,level,as.integer(p), mean))
}


