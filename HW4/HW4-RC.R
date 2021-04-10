
#function to get phat estimate

phat_f <- function(n=30,m=10000,h0=1,dist){
  #n - number of draws from distribution
  #m - number of samples drawn from distribution
  #h0 - null hypothesis
  #dist - Distribution in question, if no spec given do all distributions
  
  p<-numeric(m)
  
  #loop through the columns (individual samples) to get the estimate for that sample
  for (j in 1:m){
    
    if(dist=="norm"){
      data <- rnorm(n,h0,sd = sqrt(2))
    }else if (dist=="chi"){
      data <- rchisq(n,1)
    }else if (dist=="u"){
      data <- runif(m,0,2)
    }else if (dist=="exp"){
      data <- rexp(m,1)
    }
    ttest <- t.test(data,alternative = "two.sided",mu = h0)
    p[j] <- ttest$p.value
  }

  #this will sum up the True/(Total)
  phat.05 <- mean(p<.05)
  phat.1 <- mean(p<.1)
  #add relative error and size columns
  #keep this seperate b/c return is truncating the numbers 
  return(c("alpha05" = phat.05,"alpha1" = phat.1))
}


#distributions 
#Normal(1,2)
norm.30 <- phat_f(30,m=10000,1,dist="norm")
norm.300 <- phat_f(300,m=10000,1,dist="norm")
#Chi-Sq(1)
chi.30 <- phat_f(30,m=10000,1,dist="chi")
chi.300 <- phat_f(300,m=10000,1,dist="chi")

#Unif(0,2)
unif.30 <- phat_f(30,m=10000,1,dist="u")
unif.300 <- phat_f(300,m=10000,1,dist="u")

#exp(1)
exp.30 <- phat_f(30,m=10000,1,dist="exp")
exp.300 <- phat_f(300,m=10000,1,dist="exp")

#put all values into a table
empirical_alpha <- rbind(norm.30,norm.300,chi.30,chi.300,unif.30,unif.300,exp.30,exp.300)
