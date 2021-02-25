# Ricardo Cortez & Ben Graf
# STA 6133 - Simulation and Statistical Computing
# Homework 2
# Due 26 Feb 2021

pacman::p_load(pacman, statmod, lme4, stats)

#### 1 ####

# (a)

n <- 7
a <- 0
b <- 10
xa <- seq(from = a, to = b, by = (b-a)/(n-1))

# Set up X matrix and y vector
Xa_matrix <- matrix(nrow = n, ncol = n)
ya <- vector(mode = "numeric", length = n)
for (i in 1:n) {
  Xa_matrix[i,] <- xa^(i-1)
  ya[i] <- (b-a)^i/i
}

(wa <- solve(Xa_matrix) %*% ya)   # Find weights by solving X-inverse times y

prob1afunc <- function(x) { return( exp(-1/(1+x^2)) ) }   # Define function to be integrated

(Q1a <- sum(wa * prob1afunc(xa)))

# (b)

# Find roots of Legendre polynomial of degree 7:  1/16 * (429*x^7 - 693*x^5 + 315*x^3 - 35*x)
(xb <- sort(Re(polyroot(1/16 * c(0, -35, 0, 315, 0, -693, 0, 429)))))
Xb_matrix <- matrix(nrow = n, ncol = n)
yb <- vector(mode = "numeric", length = n)
for (i in 1:n) {
  Xb_matrix[i,] <- xb^(i-1)
  yb[i] <- (1^i-(-1)^i)/i
}
(wb <- solve(Xb_matrix) %*% yb)   # Find weights by solving X-inverse times y

#leg <- gauss.quad(n = 7, kind = "legendre")   #Alternative way to get weights and nodes

(Q1b <- (b-a)/2 * sum(wb * prob1afunc(((b-a)*xb + a + b)/2)))

# (c)

(Q1c <- integrate(prob1afunc, lower = a, upper = b))
(Q1a-Q1c$value)/Q1c$value   #Relative error between Newton-Cotes and integrate()
(Q1b-Q1c$value)/Q1c$value   #Relative error between Gauss-Legendre and integrate()

# (d)

prob1dorig <- function(x) { return( (x+1) / x^4 ) }   # Define function to be integrated
prob1dfunc <- function(t) { return( exp(t)*(t+2) / (t+1)^4 ) }   # Define function for Gauss-Laguerre

lag10 <- gauss.quad(n = 10, kind = "laguerre", alpha = 0)   #Get weights and nodes
lag100 <- gauss.quad(n = 100, kind = "laguerre", alpha = 0)   #Get weights and nodes

(Q1d_n10 <- sum(lag10$weights * prob1dfunc(lag10$nodes)))   #For n=10
(Q1d_n100 <- sum(lag100$weights * prob1dfunc(lag100$nodes)))   #For n=100

(Q1d <- integrate(prob1dorig, lower = 1, upper = Inf))
(Q1d_n10-Q1d$value)/Q1d$value   #Relative error between Gauss-Laguerre, n=10, and integrate()
(Q1d_n100-Q1d$value)/Q1d$value   #Relative error between Gauss-Laguerre, n=100, and integrate()


#### 2 ####

#### ADAPTATION OF PROFESSOR'S CODE

# We are interested in approximating the log-likelihood using numerical integration methods

# Function to compute f(beta, yi, trti, txi, sqrt(2)*sigma*u)
f <- function(u, eta, yi, trti, txi) {
  # u is a transformation of gammai; gammai = sqrt(2)*sigma*u (scalar)
  # eta is a length 5 list with beta_0, beta_1,beta_2,beta_3,sigma
  # trti is treatment indicator (scalar)
  # txi is time from the start of treatment of the jth visit for the ith patient (vector of length 7 or fewer)
  # yi is the response for the ith patient at each visit (vector of same length as txi)
  
  a <- eta[1] + eta[2]*trti + eta[3]*txi + eta[4]*txi*trti + sqrt(2)*eta[5]*u
  p <- exp(a)/(1+exp(a))
  return( prod(p^yi * (1-p)^(1-yi)) ) 
}

fv <- Vectorize(f, vectorize.args = "u")   #This allows for many values of u to be submitted to f at once 
#and evaluated in vector form, rather than a loop

# Compute the points and weights of for the Gauss-Hermite quadrature
b <- gauss.quad(n = 100, kind = "hermite")

# Import data
d <- read.table("http://faculty.business.utsa.edu/vdeolive/toenail.txt",
                header=F)
names(d) <- c("patient", "response", "treatment", "time", "visit")

# Create a patient_count variable in the dataset to denote this is the nth patient, 
# rather than using their IDs, which have gaps
cur_patient <- 0
pat_count <- 0
d$patient_count <- vector(mode = "numeric", length = nrow(d))
for (index in 1:nrow(d)) {
  if (d[index,]$patient != cur_patient) {
    cur_patient <- d[index,]$patient
    pat_count <- pat_count + 1
  }
  d[index,]$patient_count <- pat_count
}

# Function to approximate the negative log-likelihood, using Gauss-Hermite quadrature
loglik_appx <- function(eta, data) {
  # eta = (beta_0, beta_1, beta_2, beta_3, sigma)
  # data = full data set (includes at least response, treatment, time, patient_count)
  # b = dataframe with weights and nodes (Hermite), exists outside this function
  
  # Approximate the integral for each patient
  n <- length(unique(data$patient_count))
  integral_appx <- c()
  for (pat in 1:n) {
    the_rows <- which(data$patient_count == pat)
    yi <- data[the_rows,]$response
    trti <- data[the_rows[1],]$treatment
    txi <- data[the_rows,]$time
    integral_appx[pat] <- sum(b$weights * fv(b$nodes, eta, yi, trti, txi))
  }
  
  # Combine all patients with negative sum of logs
  return(-sum(log(integral_appx)))
}

# Now optimize to find the MLE
eta0 <- c(0, 0, 0, 0, 1) # Initial values
ml <- optim(eta0, loglik_appx, data = d, method = "L-BFGS-B",
            lower = c(-Inf, -Inf, -Inf, -Inf, 0), 
            upper = c(Inf, Inf, Inf, Inf, Inf), hessian = TRUE)

## Maximum likelihood estimates of beta_0, beta_1, beta_2, beta_3, and sigma
ml$par

## Approximate standard errors of the ML estimators
sqrt(diag(solve(ml$hessian)))

#### ATTEMPT TO SOLVE USING GLMER()

#Import and set up data
setwd("/Users/Ben/Library/Mobile Documents/com~apple~CloudDocs/Documents/UTSA Master's/Semester 4/STA 6133 Simulation & Statistical Computing/Homework/Homework 2")
toenail <- read.csv("./toenail.txt", header = FALSE, sep = " ")
names(toenail) <- c("patient", "response", "treatment", "time")

model <- glmer(response ~ treatment*time + (1 | patient), data = toenail, family = "binomial")
summary(model)
#plot(model)

preds <- predict(model, newdata = toenail)
binpreds <- preds>=0.5   #Bin prediction probabilities as either 0 or 1
table(toenail$response, binpreds)   #Confusion matrix
mean(toenail$response == binpreds)   #Accuracy


#### 3 ####

# (a)

my_rand <- function(m, F, interval = c(-100,100)) {
  # The function arguments must be m = the number of required draws,
  # and F, a continuous and strictly increasing cdf.
  # Added a third optional variable, interval, which is the interval over which to evaluate F (needed for root-finding).
  #   interval should be formatted as a length-2 vector as seen in the default value in the function definition.
  
  u <- runif(m)   #Draw from Uniform(0,1)
  inverse_cdf <- c()
  for (i in 1:m) {
    #Calculate inverse CDF using root-finding for each draw from Uniform(0,1)
    inverse_cdf[i] <- uniroot(function(t){F(t) - u[i]}, interval = interval)$root
  }
  return(inverse_cdf)
}


# (b)

m <- 10000

gamma_test_1 <- function(x) { return(pgamma(x, shape = 3, scale = 1)) }   #CDF of Gamma(3,1)
set.seed(1)
output_1 <- my_rand(m = m, F = gamma_test_1, interval = c(0,100))
hist(output_1, probability = TRUE, ylim = c(0,0.3), main = "Random sampling from Gamma(3,1)")
range_1 <- seq(from = 0, to = 12, by = 0.01)
lines(x = range_1, y = dgamma(range_1, shape = 3, scale = 1), col = 2)   #Add actual density curve

gamma_test_2 <- function(x) { return(pgamma(x, shape = 0.3, scale = 1)) }   #CDF of Gamma(0.3,1)
set.seed(1)
output_2 <- my_rand(m = m, F = gamma_test_2, interval = c(0,100))
hist(output_2, probability = TRUE, main = "Random sampling from Gamma(0.3,1)")
range_2 <- seq(from = 0, to = 6, by = 0.01)
lines(x = range_2, y = dgamma(range_2, shape = 0.3, scale = 1), col = 2)   #Add actual density curve

# (c)

q <- qgamma(ppoints(m), shape = 3, scale = 1)   # Quantile points for true density
qqplot(q, output_1, cex = 0.5, xlab = "Gamma(3,1)", ylab = "Sample")
abline(0,1)


#### 4 ####

# (c)

rtriangular <- function(m, a, b) {
  # The function arguments must be m = the number of required draws and the parameters a and b.
  u <- runif(m)
  Q <- ifelse(u < 0.5, 2*a + (b-a)*sqrt(2*u), 2*b - (b-a)*sqrt(2-2*u))
  return(Q)
}

a <- 2/2
b <- 8/2
set.seed(1)
samp <- rtriangular(10000, a, b)   # Simulate one random sample of size 10000 from the triangular distribution on (2, 8) = (2a, 2b)
hist(samp, probability = TRUE, ylim = c(0,0.35), main = "10,000 draws from Triangular distribution on (2,8)", xlab = "")
lines(density(samp), col = 2)

tri_pdf <- function(x, a, b) {
  # The PDF of the triangular distribution
  f <- ifelse(x<2*a | x>2*b, 0, ifelse(x<a+b, (x-2*a)/(b-a)^2, (2*b-x)/(b-a)^2) )
  return(f)
}

range <- seq(from = 2*a, to = 2*b, by = 0.01)
lines(x = range, y = tri_pdf(range, a, b), col = 4)   #Add actual PDF curve
legend("topright", legend = c("Sample density", "Actual PDF"), lty = 1, col = c(2,4), )



