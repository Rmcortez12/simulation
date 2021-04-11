library(nimble)

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


#######2
#a
# our estimators are xbar, Median, and alphaTrimmed

#need a list of theta's
#for each theta we're going to need an JXM matrix of samples from each distribution 
#J is the number of draws in the sample
#M is the simulation size

# MSE = 1/M* Sum(T(x) -g(theta))^2

#n=30, m= 100000

AllMSE.f <- function(n,theta,M){
  nM <- n*M
  
  sim_norm <- matrix(rnorm(nM,theta,9), ncol = M)
  sim_dexp <- matrix(rdexp(nM,theta,sqrt(3/2)), ncol = M)
  sim_t <- matrix(rt(nM,df=3,ncp=theta), ncol = M)  
  
  n.estimates <- matrix(estimator.f(sim_norm),ncol = 3)
  dexp.estimates <- matrix(estimator.f(sim_dexp),ncol = 3)
  t.estimates <- matrix(estimator.f(sim_t),ncol=3)
  
  #the 2 tells apply to apply the function to the columns (independent samples)
  n.mse<- apply(n.estimates,2,MSE.f,theta=theta)
  dexp.mse<-apply(dexp.estimates,2,MSE.f,theta=theta)
  t.mse <-apply(t.estimates,2,MSE.f,theta=theta)
  res <- rbind(n.mse,dexp.mse,t.mse)
  return(res)
  
}

estimator.f <- function(x){
  estimateXbar <- apply(x,2,mean)
  estimateMedian <- apply(x,2,median)
  estimateAlpha <- apply(x,2,mean,trim=0.1)
  return(c("EstimateXbar" = estimateXbar, "EstimateMedian" = estimateMedian,"EstimateAlpha" = estimateAlpha))
}
MSE.f <- function(estimate,theta){
  MSE<- mean(estimate^2)-2*theta*mean(estimate)+theta^2
  return(MSE)
}

MSEs <- AllMSE.f(30,theta = 0,100000)
colnames(MSEs) <- c("Xbar","Median","TrimmedMean")
MSEs


###2b
# looks similar but the only difference is the how the random matrix is built
# 100*(1-e) % comes from N and the rest comes from t3
n <- 30
M <- 1000

e <- .1

theta <- 0

mix.f <- function(n,M,e,theta){
  ne <- (1-e)
  a<-rnorm(n*M*ne,theta,9)
  b<-rt(n*M*e,df=3,ncp=theta)
  sim <- matrix(append(a,b),ncol=M)
  estimates <- matrix(estimator.f(sim), ncol=3)
  sim.mse <- apply(estimates,2,MSE.f,theta=theta)
return(sim.mse)
}

MixMSE.1 <- mix.f(30,100000,.1,theta)
MixMSE.3 <- mix.f(30,100000,.3,theta)
MixMSE <- rbind(MixMSE.1,MixMSE.3)
colnames(MixMSE) <- c("Xbar","Median","TrimmedMean")
MixMSE


