# Oliver Orejola
#
# Created: 7/16/21
# Free Multiplicative Convolution 
#
#
#################################

# Spectral Density of AR(1)
specdenAR <- function(v,phi,variance){
  f <- variance/(2*pi) * 1/(1-2*phi*cos(v)+phi^2)
  return(f)
}


N = 1000
U = runif(N,-pi,pi)
f <- 2*pi*diag(specdenAR(U,0.5,1))
X<- matrix(rnorm(N*N),N,N)
MP <- 1/N*crossprod(X)

E <-eigen(f %*% MP)$values
m <- max(E)

hist(E, probability=TRUE, breaks = seq(0,m+1,by = 0.05)
     #,xlim = c(0,5)
     ,xlab="Eigenvalue"
     ,main = "Hist of ESD Free Multiplicative Convolution Spectral Density of AR(1)"
)
lambda = 1
lambda_plus = (1 + sqrt(lambda))**2
lambda_minus = (1 - sqrt(lambda))**2
u=seq(lambda_minus,lambda_plus,by=.01)
v=sqrt((lambda_plus - u)*(u - lambda_minus))*1/(2*pi * u * lambda)
lines(u,v,lwd=2)

