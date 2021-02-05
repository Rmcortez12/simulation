# Ben Graf & Ricardo Cortez
# STA 6113 - Simulation
# Homework 1
# Due 9 Feb 2021

pacman::p_load(pacman, plotly)

#### 1 ####

sinapprox <- function(x, n) {
  taylor <- 0
  for (i in 1:n) {
    taylor <- taylor + (-1)^(i+1) * x^(2*i-1) / factorial(2*i-1)
  }
  return(taylor)
}
sin(pi)
sin(5*pi)
sinapprox(pi,10)
sinapprox(pi,20)
sinapprox(pi,100)
sinapprox(5*pi,10)
sinapprox(5*pi,20)
sinapprox(5*pi,100)


#### 2 ####

n <- c(1, 5, 10)
b <- 1
a <- c <- 10^(-n)
(x1 <- (-b + sqrt(b^2 - 4*a*c)) / (2*a))

polyroot(c(c[1],b,a[1]))
polyroot(c(c[2],b,a[2]))
polyroot(c(c[3],b,a[3]))

# x2 through x8 do not resolve the issue with n=10
(x2 <- (-b + sqrt((b+2*a) * (b-2*a))) / (2*a))
(x3 <- (-b + sqrt(b+2*a) * sqrt(b-2*a)) / (2*a))
(x4 <- (-b + exp((1/2)*(log(b+2*a) + log(b-2*a)))) / (2*a))
(x5 <- (-b/(2*a)) + sqrt((b^2)/(4*a^2) - c/a))
(x6 <- (-b/(2*a)) + sqrt(b/(2*a) + sqrt(c/a)) * sqrt(b/(2*a) - sqrt(c/a)))
(x7 <- (-b/sqrt(b-2*a) + sqrt(b+2*a)) * (sqrt(b-2*a) / (2*a)))
(x8 <- (-b/(2*a)) + exp((1/2)*(log(b+2*a) + log(b-2*a)) - log(2*a)))

# x9 does resolve the issue!
(x9 <- -2*c / (sqrt(b^2 - 4*a*c) + b))


#### 3(a) ####

data <- c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44, 3.29,
          3.71, -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75, 0.27, 43.21)
n <- length(data)

# Find the log-likelihood and first derivative for a range of thetas for plotting
theta <- seq(from = -60, to = 60, by = 0.1)
cauchyll <- vector(length = length(theta))
cauchydll <- vector(length = length(theta))
for (j in 1:length(theta)) {   
  cauchyll[j] <- -n*log(pi) - sum(log(1+(data-theta[j])^2))  
  cauchydll[j] <- sum(2*(data-theta[j]) / (1+(data-theta[j])^2))
}
plot(theta, cauchyll, type = "l")
plot(theta, cauchydll, type = "l")
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
    cat("At iteration", iter, "value of x is:", x, "\n")
  }
  # output depends on success of algorithm
  if (abs(fx[1]) > tol) {
    cat("Algorithm failed to converge\n")
    return(NULL)
  } else {
    cat("Algorithm converged\n")
    return(c(root = x, f.val = fx[1]))
  }
}

cauchy <- function(x) {
  # fx is first derivative of log-likelihood, the one we want the roots of
  # dfx is second derivative of log-likelihood, the first derivative of fx
  fx <- sum(2*(data-x) / (1+(data-x)^2))
  dfx <- sum(2*(-1+(data-x)^2) / (1+(data-x)^2)^2)
  return(c(fx, dfx))
}

# Test a variety of starting values
newtonraphson(cauchy, -11)
newtonraphson(cauchy, -1)
newtonraphson(cauchy, 0)
newtonraphson(cauchy, 1.5)
newtonraphson(cauchy, 4)
newtonraphson(cauchy, 4.7)
newtonraphson(cauchy, 7)
newtonraphson(cauchy, 38)
newtonraphson(cauchy, mean(data))   #This uses the sample mean as a starting point

# Find which of the Newton candidates for theta yields the highest value in the log-likelihood equation (global maximum)
-n*log(pi) - sum(log(1+(data-(-53817114377))^2))  
-n*log(pi) - sum(log(1+(data-(-0.1922866))^2))   #This is the global maximum
-n*log(pi) - sum(log(1+(data-(1.713587 ))^2))  
-n*log(pi) - sum(log(1+(data-(2.817472))^2))  
-n*log(pi) - sum(log(1+(data-(41.04085))^2))  
-n*log(pi) - sum(log(1+(data-(42.79538))^2))  
-n*log(pi) - sum(log(1+(data-(54.87662))^2))  


#### 3(b) ####

# Bisection Method function provided by Dr. DeOliveira
bisection <- function(f, x1, x2, maxit = 1000, tol = 1e-7, stop.fval = TRUE) {
  # it uses absolute error as stopping criterion when stop.fval = FALSE
  f1 <- f(x1)
  if (abs(f1) < tol)
    return(x1)
  f2 <- f(x2)
  if (abs(f2) < tol)
    return(x2)
  if (f1 * f2 > 0)
    stop("f has equal sign at endpoints of initial interval")
  if (f1 > 0) { # swap x1 and x2
    tmp <- x1 ; x1 <- x2 ; x2 <- tmp 
  }
  n <- 0 # counter
  x <- x1
  repeat {
    n <- n + 1
    x.mid <- (x1 + x2) / 2
    f.mid <- f(x.mid)
    if(stop.fval == TRUE) {
      if(abs(f.mid) < tol | n == maxit)
        break }
    else if(abs(x.mid - x) < tol | n == maxit)
      break
    if(f.mid<0){ x1 <- x.mid ; x <- x1}
    else { x2 <- x.mid ; x <- x2 }
  }
  return(list(root = x.mid, f = f.mid, iter = n))
}

cauchy2 <- function(x) {
  # first derivative of log-likelihood, which we want the roots of
  return(sum(2*(data-x) / (1+(data-x)^2)))
}

bisection(cauchy2, -1, 1, tol = 1e-9)
bisection(cauchy2, 1, 2, tol = 1e-9)   #This bracket converges to a local maximum, not the global


#### 3(c) ####

fixedpoint <- function(f, starting, alpha, maxit = 1000, tol = 1e-9) {
  # g(x) = x + alpha * f(x)
  x <- starting
  fx <- f(starting)
  iter <- 0
  while ((abs(fx) > tol) & (iter < maxit))  {
    x <- x + alpha*fx
    fx <- f(x)
    iter <- iter + 1
    cat("At iteration", iter, "value of x is:", x, "\n")
  }
  if (abs(fx) > tol) {
    cat("Algorithm failed to converge\n")
    return(NULL)
  } else {
    cat("Algorithm converged\n")
    return(list(root = x, f.val = fx, iter = iter))
  }
}

fixedpoint(cauchy2, -1, 1)
fixedpoint(cauchy2, -1, 0.64)
fixedpoint(cauchy2, -1, 0.25)


#### 4(a) ####

setwd("/Users/Ben/Library/Mobile Documents/com~apple~CloudDocs/Documents/UTSA Master's/Semester 4/STA 6133 Simulation & Statistical Computing/Homework/Homework 1/")
readLines("oilspills.dat")
oil <- read.table("oilspills.dat", header = TRUE)

# Log-likelihood function for Poisson
llpois <- function(alphas, bs, N) {
  ## alphas should be 1x2, bs should be 2xn, N should be 1xn
  term <- alphas %*% bs
  return(sum(-term + N*log(term) - log(factorial(N))))
}

# Test llpois
testalphas <- t(as.matrix(c(1,1)))
testbs <- t(as.matrix(oil[1:5,3:4]))
testN <- oil[1:5,2]
llpois(testalphas, testbs, testN)

# Gradient of log-likelihood for Poisson
gradpois <- function(alphas, bs, N) {
  ## alphas should be 1x2, bs should be 2xn, N should be 1xn
  term <- alphas %*% bs
  drv1 <- sum(-bs[1,] + N*bs[1,]/term)
  drv2 <- sum(-bs[2,] + N*bs[2,]/term)
  return(as.matrix(c(drv1,drv2)))
}

# Test gradpois
gradpois(testalphas, testbs, testN)

# Hessian of log-likelihood for Poisson
hesspois <- function(alphas, bs, N) {
  ## alphas should be 1x2, bs should be 2xn, N should be 1xn
  term <- alphas %*% bs
  drv11 <- sum(-N * (bs[1,])^2 / (term^2))
  drv12 <- sum(-N * bs[1,] * bs[2,] / (term^2))
  drv22 <- sum(-N * (bs[2,])^2 / (term^2))
  return(matrix(c(drv11,drv12, drv12,drv22), nrow = 2, ncol = 2, byrow = TRUE))
}

# Test hesspois
hesspois(testalphas, testbs, testN)

# Newton's Method adapted for this problem
newtonmulti <- function(grad, hess, x0, data1, data2, tol = 1e-9, max.iter = 100) {
  ## Newton_Raphson algorithm for solving grad == 0.
  ## x0 is the initial guess of the roots.
  ## The algorithm terminates when the function value is within distance tol of 0,
  ## or the number of iterations exceeds max.iter, whichever happens first.
  
  # initialise
  x <- x0
  fx <- grad(x, data1, data2)
  iter <-  0
  # continue iterating until stopping conditions are met
  while (((abs(fx[1]) > tol) | (abs(fx[2]) > tol)) & (iter < max.iter)) {
    x <- x - t(solve(hess(x,data1,data2)) %*% grad(x,data1,data2))
    fx <- grad(x, data1, data2)
    iter <- iter+1
    cat("At iteration", iter, "value of x is:", x, "and f is:", fx, "\n")
  }
  
  # output depends on success of algorithm
  if ((abs(fx[1]) > tol) | (abs(fx[2]) > tol)) {
    cat("Algorithm failed to converge\n")
    return(NULL)
  } else {
    cat("Algorithm converged\n")
    return(c(root = x, f.val = fx))
  }
}

# Run Newton's Method
testalphas <- t(as.matrix(c(1,1)))
bs <- t(as.matrix(oil[,3:4]))
N <- oil[,2]
res <- newtonmulti(gradpois, hesspois, testalphas, bs, N)
llpois(res[1:2], bs, N)

newtonmulti(gradpois, hesspois, t(as.matrix(c(0.5,0.5))), bs, N)
newtonmulti(gradpois, hesspois, t(as.matrix(c(0.1,0.5))), bs, N)
newtonmulti(gradpois, hesspois, t(as.matrix(c(0.8,0.5))), bs, N)
newtonmulti(gradpois, hesspois, t(as.matrix(c(3,0.5))), bs, N)
newtonmulti(gradpois, hesspois, t(as.matrix(c(2,0.5))), bs, N)
newtonmulti(gradpois, hesspois, t(as.matrix(c(2,2))), bs, N)


# Plot the log-likelihood for a variety of alphas to check that we are converging to a maximum
# Try alphas from 0 to 10
a1 <- a2 <- seq(from = 0, to = 10, by = 0.1)
len <- length(a1)^2
llik <- data.frame(alpha1 = rep_len(0,length.out = len), alpha2 = rep_len(0,length.out = len), value = rep_len(0,length.out = len))
index <- 0
# Plot log-likelihood
for (i in 1:length(a1)) {
  for (j in 1:length(a2)) {
    index <- index + 1
    val <- llpois(c(a1[i],a2[j]), bs, N)
    llik[index,] <- c(a1[i], a2[j], val)
  }
}
plot_ly(llik, x = ~alpha1, y = ~alpha2, z = ~value)
max(llik$value, na.rm = TRUE)
llik[which(llik$value == max(llik$value, na.rm = TRUE)),]


#### 4(b) ####

# Fisher Information of log-likelihood for Poisson
ipois <- function(alphas, bs, N) {
  ## alphas should be 1x2, bs should be 2xn, N should be 1xn
  term <- alphas %*% bs
  drv11 <- -sum((bs[1,])^2 / term)
  drv12 <- -sum(bs[1,] * bs[2,] / term)
  drv22 <- -sum((bs[2,])^2 / term)
  return(matrix(c(drv11,drv12, drv12,drv22), nrow = 2, ncol = 2, byrow = TRUE))
}

# Test ipois
ipois(testalphas, testbs, testN)

# Run Newton's Method
res2 <- newtonmulti(gradpois, ipois, testalphas, bs, N)
llpois(res2[1:2], bs, N)

newtonmulti(gradpois, ipois, t(as.matrix(c(0.5,0.5))), bs, N)
newtonmulti(gradpois, ipois, t(as.matrix(c(0.1,0.5))), bs, N)
newtonmulti(gradpois, ipois, t(as.matrix(c(0.8,0.5))), bs, N)
newtonmulti(gradpois, ipois, t(as.matrix(c(3,0.5))), bs, N)
newtonmulti(gradpois, ipois, t(as.matrix(c(2,0.5))), bs, N)
newtonmulti(gradpois, ipois, t(as.matrix(c(2,2))), bs, N)


#### 4(c) ####

# Run time to generate Hessian matrix
start_time <- Sys.time()
hesspois(testalphas, testbs, testN)
end_time <- Sys.time()
print(paste0("Hessian method takes ", end_time-start_time))

# Run time to generate Fisher Information matrix
start_time <- Sys.time()
ipois(testalphas, testbs, testN)
end_time <- Sys.time()
print(paste0("Fisher Information method takes ", end_time-start_time))


#### 4(d) ####

sqrt(diag(solve(-hesspois(res[1:2], bs, N))))   # Estimated standard errors of MLEs is square root of diagonal of inverse Hessian matrix
sqrt(diag(solve(-ipois(res[1:2], bs, N))))   # Estimated standard errors of MLEs is square root of diagonal of inverse Fisher Information matrix
