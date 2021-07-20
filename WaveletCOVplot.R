# Oliver Orejola
# 6/30/21
#
# Wavelet PLAY

library(somebm)
library(wavelets)

H = 0.7
j=20
X<- fbm(hurst = H ,2^j)
wave <- dwt(X, filter ="d4", n.levels = j)
WW <- wave@W
WV <- wave@V
XX <- c()
for(i in 1:wave@level){
  XX <- append(XX, log(sum(unlist(WW[i])^2)/length(WW[i])))
}
plot(XX, ylab = "log Wavelet Sample Cov", xlab = "scale")
#u=seq(0,j,by=.01)
#v=(2*H+1)*u*log(2)
#lines(u,v,lwd=2)
