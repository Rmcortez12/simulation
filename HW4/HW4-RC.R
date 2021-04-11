####1a
#function to get phat estimate

phat_f <- function(n=30,m=10000,h0=1,dist, mu1 =1){
  #n - number of draws from distribution
  #m - number of samples drawn from distribution
  #h0 - null hypothesis
  #dist - Distribution in question, if no spec given do all distributions
  
  p<-numeric(m)
  
  #loop through the columns (individual samples) to get the estimate for that sample
  for (j in 1:m){
    
    if(dist=="norm"){
      data <- rnorm(n,mu1,sd = sqrt(2))
    }else if (dist=="chi"){
      data <- rchisq(n,mu1)
    }else if (dist=="u"){
      data <- runif(m,0,2*mu1)
    }else if (dist=="exp"){
      data <- rexp(m,mu1)
    }
    ttest <- t.test(data,alternative = "two.sided",mu = h0)
    p[j] <- ttest$p.value
  }

  #this will sum up the True/(Total)
  phat.05 <- mean(p<.05)
  phat.1 <- mean(p<.1)
  #add relative error and size columns
  #keep this seperate b/c return is truncating the numbers 
  
  #power for number 1b
  
  
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


#####1B
#Power is the probability of rejection given a mu
#  p<.05 
# re-use above function, modify to take in different mu

mus <- c(seq(0,4,.25))
l <- length(mus)
n.power<-numeric(l)
c.power<-numeric(l)
u.power<-numeric(l)
e.power<-numeric(l)
for(i in 1:l){
  n.power[i] <- phat_f(300,m=10000,1,dist="norm",mu1 = mus[i])["alpha05"]
  c.power[i] <- phat_f(300,m=10000,1,dist="chi",mu1 = mus[i])["alpha05"]
  u.power[i] <- phat_f(300,m=10000,1,dist="u",mu1 = mus[i])["alpha05"]
  #e.power[i] <- phat_f(300,m=10000,1,dist="exp",mu1 = mus[i])["alpha05"]
}
for(i in 2:l){
  e.power[i] <- phat_f(300,m=10000,1,dist="exp",mu1 = mus[i])["alpha05"]
}
plot(mus,n.power)
plot(mus,c.power)
plot(mus,u.power)
plot(mus,e.power)
