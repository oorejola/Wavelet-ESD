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
#################################

N = 750
U = runif(N,-pi,pi)
phi = 0.5
f <- 2*pi*diag(specdenAR(U,phi,1))
X<- matrix(rnorm(N*N),N,N)
MP <- 1/N*crossprod(X)

X<-c()
for(i in 1:N){
  X<- append(X,arima.sim(list(order=c(1,0,0), ar=phi), n=N))
}
XX <- matrix( X, N,N, byrow = TRUE)  
E_AR <-eigen(1/N* crossprod(XX))$values

E_FMP <-eigen(f %*% MP)$values
M <- max(E_AR,E_FMP)


par(mfrow = c(2,1))
hist(E_FMP, probability=TRUE, breaks = seq(0,M+1,by = 0.05)
     ,ylim=c(0,2)
     ,xlab= "Eigenvalues"
     ,main = paste("ESD FMC Spectral Density of AR(1) phi =", phi)
)
addMP(1)
hist(E_AR,probability=TRUE
     ,breaks = seq(0,M+1, by = 0.05)
     ,ylim=c(0,2)
     ,xlab= "Eigenvalues"
     ,main = paste("ESD Sample Cov AR(1) phi =", phi)
)
addMP(1)


