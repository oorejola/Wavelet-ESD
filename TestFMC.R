# Oliver Orejola
#
# Created: 7/16/21
# Updated: 7/27/21
# Free Multiplicative Convolution 
#
#
#################################

addMP<- function(lambda){
  lambda_plus = (1 + sqrt(lambda))**2
  lambda_minus = (1 - sqrt(lambda))**2
  u=seq(lambda_minus,lambda_plus,by=.01)
  v=sqrt((lambda_plus - u)*(u - lambda_minus))*1/(2*pi * u * lambda)
  lines(u,v,lwd=2)
}

# Spectral Density of AR(1)
specdenAR <- function(v,phi,variance){
  f <- variance/(2*pi) * 1/(1-2*phi*cos(v)+phi^2)
  return(f)
}
# Spectral Density of fARIMA(0,d,0)
specdenfARIMA <- function(v,d){
  f<- (1/(2*pi))*((2*sin(v/2))^(2))^(-d)
  return(f)
}
#################################

N = 750
U = runif(N,-pi,pi)
phi = 0.75
d = -0.2
f_AR <- 2*pi*diag(specdenAR(U,phi,1))

f_fARIMA <- 2*pi*diag(specdenfARIMA(U,d))
X<- matrix(rnorm(N*N),N,N)
MP <- 1/N*crossprod(X)

E_AR_FMP <-eigen(f_AR %*% MP)$values
E_fARMA_FMP <-eigen(f_fARIMA %*% MP)$values

##############################
X<-c()
for(i in 1:N){
  X<- append(X,arima.sim(list(order=c(1,0,0), ar=phi), n=N))
}
XX <- matrix( X, N,N, byrow = TRUE)  
E_AR <-eigen(1/N* crossprod(t(XX)))$values

M <- max(E_AR,E_AR_FMP)

par(mfrow = c(2,1))
hist(E_AR_FMP, probability=TRUE, breaks = seq(0,M+1,by = 0.05)
     ,ylim=c(0,2)
     ,xlab= "Eigenvalues"
     ,main = paste("ESD via free multiplicative convolution  AR(1) phi=", phi)
)
addMP(1)
hist(E_AR,probability=TRUE
     ,breaks = seq(0,M+1, by = 0.05)
     ,ylim=c(0,2)
     ,xlab= "Eigenvalues"
     ,main = paste("ESD Sample Cov AR(1) phi =", phi)
)
addMP(1)
#################################

X<- c()
for(i in 1:N){
  fARMA<- simARMA0(N, H = d+0.5)
  X <- append(X,fARMA)
}
XX <- matrix( X, N, N,byrow = TRUE)  
E_fARMA <-eigen(1/N*crossprod(t(XX)))$values
M <- max(E_fARMA,E_fARMA_FMP)

hist(E_fARMA_FMP, probability=TRUE, breaks = seq(0,M+1,by = 0.05)
     ,ylim=c(0,2)
     ,xlab= "Eigenvalues"
     ,main = paste("ESD via free multiplicative convolution  fARIMA d =", d)
)
addMP(1)
hist(E_fARMA, probability=TRUE, breaks = seq(0,M+1,by = 0.05)
     ,ylim=c(0,2)
     ,xlab= "Eigenvalues"
     ,main = paste("ESD Sample Cov fARIMA d =", d)
)
addMP(1)