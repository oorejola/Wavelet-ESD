# Oliver Orejola
# Created: 6/30/21
# Updated: 7/17/21
# Wavelet ESD

#########[ Run First ]#########
library(somebm)
library(wavelets)

waveESD <-function(level,scale,H){
  N <- 2^level
  n_j <- N/(2^scale)
  #factor <- (2*H+1)*scale*log(2)
  X <- c()
  for(i in 1:n_j){
    fbmSIM <- N^H* fbm(hurst = H, n = N)
    wave <- dwt(fbmSIM, filter = "d4", n.levels = level, boundary = "reflection")
    X<- append(X,unlist(wave@W[scale]))
  }
  XX <- matrix( X, n_j, n_j,byrow = TRUE)  
  E <-log(eigen(1/n_j* crossprod(XX))$values)#-factor
  return(E)
}

addlogMP <-function(lambda){
  lambda_plus = (1 + sqrt(lambda))**2
  lambda_minus = (1 - sqrt(lambda))**2
  u=seq(-50,log(4),by=.0005)
  v=sqrt((lambda_plus - exp(u))*( exp(u) - lambda_minus))*1/(2*pi  * lambda)
  lines(u,v,lwd=2)
}
###############################
level <- 12
scale <- 6
Hurst <- 0.1

E<- waveESD(level,scale,Hurst)
m <- min(E)
M <- max(E)

hist(E)
hist(E,probability=TRUE
     ,breaks = seq(m-1,M+1, by = 0.3)
     ,xlim=c(m,M+1) ,ylim=c(0,0.4), xlab= "Eigenvalues"
     ,main = paste("H=",Hurst,", n=2^",level,", j=",scale)
)
addMP(1)

##############################
par(mfrow = c(5,4))
level <- 20

H_list <- c(0.1,0.25,0.5,0.75,0.9)
scales <- c(12,14,16,18)

for(H in H_list){
  E1 <- waveESD(level, scale = scales[1], H)
  E2 <- waveESD(level, scale = scales[2], H)
  E3 <- waveESD(level, scale = scales[3], H)
  E4 <- waveESD(level, scale = scales[4], H)
  
  m <- min(E1,E2,E3,E4)
  M <- max(E1,E2,E3,E4)

  hist(E1,probability=TRUE
       ,breaks = seq(m-1,M+1, by = 0.3)
       ,xlim=c(m,M) ,ylim=c(0,0.4), xlab= "Eigenvalues"
       ,main = paste("H=",H,", n=2^",level,", j=",scales[1])
      )
  addlogMP(1)

  hist(E2,probability=TRUE
       ,breaks = seq(m-1,M+1, by = 0.5)
       ,xlim=c(m,M) ,ylim=c(0,0.4), xlab= "Eigenvalues"
       ,main = paste("H=",H,", n=2^",level,", j=",scales[2])
  )
  addlogMP(1)
  
  hist(E3,probability=TRUE
       ,breaks = seq(m-1,M+5, by = 0.8)
       ,xlim=c(m,M) ,ylim=c(0,0.4), xlab= "Eigenvalues"
       ,main = paste("H=",H,", n=2^",level,", j=",scales[3])
  )
  addlogMP(1)
  hist(E4,probability=TRUE
     ,breaks = seq(m-1,M+5, by = 1.5)
     ,xlim=c(m,M) ,ylim=c(0,0.4), xlab= "Eigenvalues"
     ,main = paste("H=",H,", n=2^",level,", j=",scales[4])
  )
  addlogMP(1)
}

###
# fixed MP
###
par(mfrow = c(5,4))
level <- 20

H_list <- c(0.1,0.25,0.5,0.75,0.9)
scales <- c(12,14,16,18)

for(H in H_list){
  E1 <- waveESD(level, scale = scales[1], H)
  E2 <- waveESD(level, scale = scales[2], H)
  E3 <- waveESD(level, scale = scales[3], H)
  E4 <- waveESD(level, scale = scales[4], H)
  
  m <- min(E1,E2,E3,E4)
  M <- max(E1,E2,E3,E4)
  
  hist(E1,probability=TRUE
       ,breaks = seq(m-1,M+1, by = 0.3)
       ,xlim=c(-15,35) ,ylim=c(0,0.4), xlab= "Eigenvalues"
       ,main = paste("H=",H,", n=2^",level,", j=",scales[1])
  )
  addlogMP(1)
  
  hist(E2,probability=TRUE
       ,breaks = seq(m-1,M+1, by = 0.5)
       ,xlim=c(-15,35) ,ylim=c(0,0.4), xlab= "Eigenvalues"
       ,main = paste("H=",H,", n=2^",level,", j=",scales[2])
  )
  addlogMP(1)
  
  hist(E3,probability=TRUE
       ,breaks = seq(m-1,M+5, by = 0.8)
       ,xlim=c(-15,35) ,ylim=c(0,0.4), xlab= "Eigenvalues"
       ,main = paste("H=",H,", n=2^",level,", j=",scales[3])
  )
  addlogMP(1)
  hist(E4,probability=TRUE
       ,breaks = seq(m-1,M+5, by = 1.5)
       ,xlim=c(-15,35) ,ylim=c(0,0.4), xlab= "Eigenvalues"
       ,main = paste("H=",H,", n=2^",level,", j=",scales[4])
  )
  addlogMP(1)
}


#H = 0.25
#k=10
#fBM<-((2^k)^H)*fbm(hurst=H,n=2^k)
#plot(fBM)
