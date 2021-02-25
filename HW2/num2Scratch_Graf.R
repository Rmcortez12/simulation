library(statmod)

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
