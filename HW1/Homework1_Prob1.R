#1 taylor expansion 

#Function for taylor expansion of sin(x)
taylor_expansion <- function(n,x){
  p <- (-1)^(n+1)*((x^(2*n-1))/factorial(2*n-1))
  return(p)
}

#Function for summation of expansion
loop_func <- function(n,x){
  sum = 0 
  for(i in 1:n){
    p = taylor_expansion(i,x)
    sum = sum+p
  }
  return(sum)
}

#Trying to consolidate data into data frame

sol <- data.frame("N"=c(10,20,100,10,20,100),"X"=c(pi,5*pi,pi,5*pi,pi,5*pi))

#x = pi (sin(PI)=0)
loop_func(10,pi)
loop_func(20,pi)
loop_func(100,pi)

#x = 5*pi (sin(5*PI)=0)
loop_func(10,5*pi)
loop_func(20,5*pi)
loop_func(100,5*pi)
loop_func(120,5*pi)


####2####
quad <- function(n,b){
  a <- 10^-n
  c <- a
  root_num_right <- sqrt(b^2-(4*a*c))
  print(root_num_right)
  root_num <- -b +root_num_right
  root_den <- 2*a
  return(root_num/root_den)
}

quadTest <- function(n,b){
  a <- 10^-n
  c <- a
  root_num = -4*a*c
  root_den = sqrt(b^2-(4*a*c))+b
  return(root_num/root_den)
}

quadTest(1,1)
quadTest(5,1)
quadTest(10,1)


quad(1,1)
quad(5,1)
quad(10,1)


####3a####
data <- c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44, 3.29,
          3.71, -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75, 0.27, 43.21)

n <- length(data)

log_cauchy <- function(theta){
  return(n*log(pi)-sum(log(1+(data-theta)^2)))
}

log_cauchy_prime <- function(theta){
return(-sum(2*(data-theta)/(1+(data-theta)^2)))
}

log_cauchy_2prime <- function(theta){
  return(-sum(2*(-1+(data-theta)^2) / (1+(data-theta)^2)^2))
}

cauchy <- function(theta){
  return(c(log_cauchy(theta),
           log_cauchy_prime(theta)))
}

cauchy_primes <-function(theta){
  return(c(log_cauchy_prime(theta),log_cauchy_2prime(theta)))
}

theta_list <- seq(-50,50,.1)
values = sapply(theta_list, cauchy)

log_cauchy_values <- values[1,]
log_cauchy_prime_values <- values[2,]

plot(theta_list,log_cauchy_values, type = "l")
plot(theta_list,log_cauchy_prime_values, type = "l")
abline(h = 0, col = "red")

# Newton's Method function provided by Dr. DeOliveira
newtonraphson <- function(ftn, x0, tol = 1e-9, max.iter = 100) {
  ## Newton_Raphson algorithm for solving ftn(x)[1] == 0
  ## It is assumed that ftn is a function of a single variable that returns
  ## the function value and its first derivative as a vector of length 2.
  ## x0 is the initial guess of the root.
  ## The algorithm terminates when the function value is within distance tol of 0,
  ## or the number of iterations exceeds max.iter, whichever happens first
  # initialise
  x <- x0
  fx <- ftn(x)
  iter <-  0
  # continue iterating until stopping conditions are met
  while ((abs(fx[1]) > tol) & (iter < max.iter)) {
    x <- x - fx[1]/fx[2]
    fx <- ftn(x)
    iter<- iter+1
    #cat("At iteration", iter, "value of x is:", x, "\n")
  }
  # output depends on success of algorithm
  if (abs(fx[1]) > tol) {
    #cat("Algorithm failed to converge\n")
    return(NULL)
  } else {
    #cat("Algorithm converged\n")
    return(c(root = x, f.val = fx[1]))
  }
}

mle_values <- roots <- test_values <- c(-11,-1,0,1.5,4,4.7,7,38)

for(i in 0:length(test_values)){
  root = newtonraphson(cauchy_primes,test_values[i])[1]
  roots[i] = root
  mle_values[i] = log_cauchy(root)
}
roots
mle_values
test_values[which.max(mle_values)] ### max value is at -1
