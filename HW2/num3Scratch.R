set.seed(10000)
x <- rnorm(1000, 10, 13)
pdf <- density(x)

# Interpolate the density
f <- approxfun(pdf$x, pdf$y, yleft=0, yright=0)
# Get the cdf by numeric integration
cdf <- function(x){
  integrate(f, -Inf, x)$value
}
# Use a root finding function to invert the cdf
invcdf <- function(q){
  uniroot(function(x){cdf(x) - q}, a)$root
}

med <- invcdf(.5)
cdf(med)



F <- function(x) 1-(1+x)*exp(-x)

F.inv <- function(y){uniroot(function(x){F(x)-y},interval=c(0,100))$root}
F.inv <- Vectorize(F.inv)

x <- seq(0,10,length.out=1000)

X <- runif(1000,0,1)   # random sample from U[0,1]
Z <- F.inv(X)


#################################
#recall cdf is just integral of pdf
#given a cdf F, set equal unif() then solve for x
#There exists a real number for which F(x) = y
# so F(x) - y = 0 would give us that real Number

sim <- function(F,m){
  x <- seq(0,m)
  
  F.inv <- function(y){uniroot(function(x){F(x)-y},interval=c(0,100))$root}
  F.inv <- Vectorize(F.inv)

  X <- runif(1000,0,1)   # random sample from U[0,1]
  Z <- F.inv(X)
  return(Z)
}
sim(F,1000)

