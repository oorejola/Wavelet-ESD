#Oliver Orejola
#Created: 7/19/21
# QQ-Plot

#########[ Run First ]#########
library(somebm)
library(wavelets)
library(RMTstat)
library(ggplot2)

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

addlogMP <-function(){
  lambda_plus = (1 + sqrt(1))**2
  lambda_minus = (1 - sqrt(1))**2
  u=seq(-50,log(4),by=.0005)
  v=sqrt((lambda_plus - exp(u))*( exp(u) - lambda_minus))*1/(2*pi )
  lines(u,v,lwd=2)
}

addMP <- function(lambda){
  lambda_plus = (1 + sqrt(lambda))**2
  lambda_minus = (1 - sqrt(lambda))**2
  u=seq(lambda_minus,lambda_plus,by=.01)
  v=sqrt((lambda_plus - u)*(u - lambda_minus))*1/(2*pi * u * lambda)
  lines(u,v,lwd=2)
}
###############################
par(mfrow = c(3,1))

level <- 11
scale <- 2
Hurst <- 0.1

E<- exp(waveESD(level,scale,Hurst))
m <- min(E)
M <- max(E)

hist(E,probability=TRUE
     ,breaks = seq(m-1,M+1, by = 0.05)
     , xlab= "Eigenvalues"
     ,xlim=c(m,M+1) #,ylim=c(0)
     ,main = paste("H=",Hurst,", n=2^",level,", j=",scale)
)
addMP(1)

addMP(1)
qmp( 5, ndf=NA, pdim=1, var=1, svr=1, lower.tail = TRUE, log.p = FALSE )
#X<-rmp( 2000, ndf=NA, pdim=1, var=1, svr=1 )
#hist(X,probability=TRUE
#     ,breaks = seq(0,5,by=0.05))



quantile(mp)

set.seed(123);
#DATA <- rexp(40,1);
DATA <- E;


#QQ plot of data against MP(1) distribution
N         <- length(DATA);
PERCS     <- ((1:N)-0.5)/N;
#QUANTILES <- -log(1-PERCS);
QUANTILES <-qmp(PERCS, ndf=NA, pdim=1, var=1, svr=1, lower.tail = TRUE, log.p = FALSE )
PLOTDATA <- data.frame(Sample = sort(DATA),
                       Theoretical = QUANTILES);


#Generate custom QQ plot

theme_update(plot.title    = element_text(size = 15, hjust = 0.5),
             plot.subtitle = element_text(size = 10, hjust = 0.5),
             axis.title.x  = element_text(size = 10, hjust = 0.5),
             axis.title.y  = element_text(size = 10, vjust = 0.5),
             plot.margin   = unit(c(1, 1, 1, 1), "cm"));

QQPLOT <- ggplot(data = PLOTDATA, aes(x = Theoretical, y = Sample)) +
  geom_point(size = 2, colour = 'blue') +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') + 
  ggtitle(paste('Quantile-Quantile Plot H=',Hurst)) + 
  labs(subtitle = '(Comparison to Marchenko Pastur)') + 
  xlab('Theoretical Quantiles') + 
  ylab('Sample Quantiles');

QQPLOT



