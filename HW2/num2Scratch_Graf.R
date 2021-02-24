
# We are intested in approximating the log-likelihood using numerical integration methods

#Function to compute f(beta, yi, sqrt(2)*sigma*u)
f <- function(u, eta, yi, trti, txi) {
  #u is 
  #eta is a length 5 list with beta_0, beta_1,beta_2,beta_3,sigma
  #trti is treatment indicator [list] (294 elements)
  #txi is time from the start of treatment of the jth visit for the ith patient is a jxi matrix (7 rows x 294 columns)

  a <- eta[1] + eta[2]*trti + eta[3]*txi + eta[4]*txi*trti + sqrt(2)*eta[5]*u
  p <- exp(a)/(1+exp(a))
  return( prod(p^yi * (1-p)^(1-yi)) ) 
}

fv <- Vectorize(f, vectorize.args = "u")

#bring in necessary libraries
library(statmod)
## compute the points and weights of for the Gauss-Hermite quadrature
##why 100??
b <- gauss.quad(n = 100, kind = "hermite")

## Function to approximate the log-likelihood l_i(eta, yi),
## up to a constant [eta], using Gauss-Hermite quadrature
## eta = (beta_0, beta_1, sigma); yi = data from patient i; b = dataframe with weights and nodes (hermite)
#remember first we approximate the integral inside of the log, and then take the log of that value

d <- read.table("http://faculty.business.utsa.edu/vdeolive/toenail.txt",
                header=F)
names(d) <- c("patient", "response", "treatment", "time", "visit")
#patient ID, binary response, treatment code, time of eval
#trt <- d[,3]


# Create a patient_count variable in the dataset to denote this is the nth patient 
# (rather than using their IDs, which have gaps)
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


# Negative log-likelihood approximation function
loglik_appx <- function(eta, data){
  
  # Create 7x294 matrix of times, 294x7 data frame of responses (y), and 294-length vector of treatments (trt)
  tx <- matrix(nrow = 7, ncol = 294)
  y <- data.frame(visit1 = integer(),
                  visit2 = integer(), 
                  visit3 = integer(), 
                  visit4 = integer(), 
                  visit5 = integer(), 
                  visit6 = integer(), 
                  visit7 = integer() )
  trt <- vector(mode = "integer", length = ncol(tx))
  for (index in 1:nrow(data)) {
    tx[ data[index,]$visit, data[index,]$patient_count ] <- data[index,]$time
    y[ data[index,]$patient_count, data[index,]$visit ] <- data[index,]$response
    if (data[index,]$visit == 1) { trt[ data[index,]$patient_count ] <- data[index,]$treatment }
  }
  
  # Now approximate the integral for each patient
  integral_appx <- vector(mode = "numeric", length = nrow(y))
  for (pat in 1:nrow(y)) {
    integral_appx[pat] <- sum(b$weights * f(b$nodes, eta, y[pat,], trt[pat], tx[,pat]))
  }
  
  # Combine all patients with negative log product
  return(-log(prod(integral_appx)))
}


loglik_appx_2 <- function(eta, data) {
  
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
  
  # Combine all patients with negative log product
  return(-sum(log(integral_appx)))
}



eta0 <- c(0, 0, 0, 0, 1) # Initial values
ml <- optim(eta0, loglik_appx_2, data = d, method = "L-BFGS-B",
            lower = c(-Inf, -Inf, -Inf, -Inf, 0), 
            upper = c(Inf, Inf, Inf, Inf, Inf), hessian = TRUE)

## Maximum likelihood estimates of beta_0, beta_1, beta_2, beta_3, and sigma
ml$par

## Approximate standard errors of the ML estimators
sqrt(diag(solve(ml$hessian)))
