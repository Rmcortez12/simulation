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


#####3
library(lawstat) #needed for levene.test
n1<-n2 <- 10
n3<-20
alpha <- 0.05
m <- 10000
#dists
#normal(0,sig)
#exp()


#bart's test is one x,y where x is all values and y is the sample they come from
library(data.table)
samp1<- data.table("x" = rnorm(10,0,1), "y" = "samp1")
samp2<- data.table("x" = rnorm(10,0,1), "y" = "samp2")
samp3<- data.table("x" = rnorm(20,0,1), "y" = "samp3")
all <- rbind(samp1,samp2,samp3)
b<- bartlett.test(x~y, data = all)
b$p.value
#this is oddly similar to number 1??


phat.f3 <- function(n,m=10000,dist,std){
  #n - number of draws from distribution [vector]
  #m - number of samples drawn from distribution
  #dist - Distribution in question, if no spec given do all distributions
  #std - null std deviation values
  #std1 - alt std deviationvalues
  
  p<-numeric(m)
  
  #loop through the columns (individual samples) to get the estimate for that sample
  for (j in 1:m){
    simdata <- data.table(x=numeric(),y=numeric())
    for(i in 1:length(n)){
      if(dist=="norm"){
        data<-rnorm(n[i],0,sqrt(std[i]))
      }else if (dist=="t"){
        data <- rt(n[i],3,0)*sqrt(std[i])
      }else if (dist == "e"){
        data <- rexp(n[i],sqrt(std[i]))
      }
      simdata <- rbind(data.table(x = data,y = i),simdata)
    }
    b <- bartlett.test(x~y, data = simdata)
    p[j] <- b$p.value
  }
  
  #this will sum up the True/(Total)
  phat.05 <- mean(p<.05)
  #add relative error and size columns
  #keep this seperate b/c return is truncating the numbers 
  
  #power for number 1b
  
  
  return(c("alpha" = phat.05))
}

n.null.1 <- phat.f3(c(10,10,20),10000,"norm",c(1,1,1))
n.null.10 <- phat.f3(c(10,10,20),10000,"norm",c(10,10,10))

t.null.1 <- phat.f3(c(10,10,20),10000,"t",c(1,1,1))
t.null.10 <- phat.f3(c(10,10,20),10000,"t",c(10,10,10))

e.null.1 <- phat.f3(c(10,10,20),10000,"e",c(1,1,1))
e.null.10 <- phat.f3(c(10,10,20),10000,"e",c(10,10,10))

(sig_levels <- rbind(n.null.1,n.null.10,t.null.1,t.null.10,e.null.1,e.null.10))


#Power

n.alt1 <- phat.f3(c(10,10,20),10000,"norm",c(1,1,3))
n.alt2 <- phat.f3(c(10,10,20),10000,"norm",c(1,2,6))

t.alt1 <- phat.f3(c(10,10,20),10000,"t",c(1,1,3))
t.alt2 <- phat.f3(c(10,10,20),10000,"t",c(1,2,6))

e.alt1 <- phat.f3(c(10,10,20),10000,"e",c(1,1,3))
e.alt2 <- phat.f3(c(10,10,20),10000,"e",c(1,2,6))

(sig_levels <- rbind(n.alt1,n.alt2,t.alt1,t.alt2,e.alt1,e.alt2))



######4
data.4 <- c(3,5,7,18,43,85,91,98,100,130,230,487)


#We are drawing these bootstrap samples from an exp(u) sample where u is the rate MLE of exp (xbar)
#   remember though we need 'B' data samples each containing n samples 
xbar <- mean(data.4)
n <- length(data.4)
B <- 10000
set.seed(1)
sampleData <- rexp(n*B,1/xbar)
#put this sample data into matrix format that breaks all the data into B samples made of n entries
sampleMat <- matrix(sampleData,nrow = n,ncol = B)

T.samp <- 1/apply(sampleMat, 2, mean)   # M draws from bootstrap distribution of g(mu)=1/mu
hist(T.samp)
(bias_est <- mean(T.samp) - (1/xbar))
(sd_est <- sd(T.samp))

#b 
# Exact confidence intervals
# xbar dist is gamma(n,u/n)
(low <- xbar/qgamma(.90,n,1/n))
(high <- xbar/qgamma(.10,n,1/n))
(a.value <- sum(sampleData>5,na.rm = T)/length(sampleData))

#Parametric Asymptotic 
(low.a <- xbar - qnorm(.1)*xbar/sqrt(n))
(high.a <- xbar + qnorm(.1)*xbar/sqrt(n))


#bootstrap
v <- quantile(T.samp, probs = c(0.9, 0.1))
c(2*mean(data.4) - v[1], 2*mean(data.4) - v[2])
