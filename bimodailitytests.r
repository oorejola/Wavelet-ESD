

library(diptest)

bimodal_dist <- function(n,mode_1,mode_2,var_1,var_2){
  data <- c()
  for(i in c(1:n)){
    data[i]<- sample(c(rnorm(1,mode_1,sqrt(var_1)),rnorm(1,mode_2,sqrt(var_2))),1)
  }
  return(data)
}

test<- bimodal_dist(2^6,0.4,0.5,0.001,0.001)
dip.test(test)$p.value

hist(test, breaks =32,xlim = c(0,1),main = "Bimodal Hist")

monte_carlo_bimodal <- function(runs,dip_test_significance,n,mode_1,mode_2,var_1,var_2){
  count <-  0
  for(i in 1:runs){
    if(i %% 100 ==0){
      print(paste("run :", i))
      print(count)
    }
    estimate <- bimodal_dist(n,mode_1,mode_2,var_1,var_2)
    pvalue<-dip.test(estimate)$p.value
    if( pvalue<= dip_test_significance){
    count <- count+ 1
     #print("reject")
    }
  }
  proportion_of_rejections <- count/runs
  return(proportion_of_rejections)
}

monte_carlo_bimodal(1000,0.65,2^6,0.4,0.5,0.1,0.1)

runs = 1000
sig_levels <-  c(seq(0.05,0.95,by = 0.05),seq(0.96,1,by = 0.01))
alphas_bimodal <-c()
index <- 1

for(sig in sig_levels){
  test <- monte_carlo_bimodal(runs,sig,2^6,0.4,0.5,0.001,0.001)
  print(paste("dip_test_significance:" ,sig))
  print(test)
  alphas_bimodal[index] <- test
  index <- index + 1 
}

plot(sig_levels,alphas_bimodal,type = "o",
     xlab="dip test significance", ylab= "MC rejection proprtion")

hist(bimodal_dist(2^16,0.4,0.5,0.001,0.001), breaks =32,xlim = c(-1,1),main = "Bimodal Hist")