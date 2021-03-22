#HW3

#1 
p <- function(n,alpha,beta){
  
  u<-runif(10000)
  
  b <- rbeta(u,alpha,beta)
  #print(b)
  c <- factorial(n)/(factorial(b)*factorial(n-b))
  print(c)
  p <- gamma(alpha+beta)*gamma(alpha+b)*gamma(n-b+beta)/(gamma(alpha)*gamma(beta)*gamma(n+alpha+beta))
  return(p*c)
}

#testing to match wikipedia <- https://en.wikipedia.org/wiki/Beta-binomial_distribution#/media/File:Beta-binomial_distribution_pmf.png
alpha <- 0.2
beta <-.25

hist(p(20,alpha,beta))


alpha <- 3
beta <-2

hist(p(20,alpha,beta))



