# Ben Graf & Ricardo Cortez
# STA 6113 - Simulation and Statistical Computing
# Homework 3
# Due 26 Mar 2021

pacman::p_load(pacman)

#### 1 ####

# (b)

# Function to make one draw from a beta-binomial distribution
rbeta_bin <- function(n, alpha, beta) {
  p <- rbeta(1, alpha, beta)
  x <- sample(c(0,1), size = n, replace = TRUE, prob = c(1-p, p))
  return(sum(x))
}

# (c)

n <- 20
alpha <- 3
beta <- 2
iter <- 10000

# Make 10000 draws from beta-binomial
samp <- c()
set.seed(1)
for (j in 1:iter) {
  samp[j] <- rbeta_bin(n, alpha, beta)
}
mean(samp)  # Sample mean
var(samp)   # Sample variance
hist(samp, probability = TRUE)   # Sample histogram
lines(density(samp))   # Sample density

# Compare with theoretical
range <- 0:20
# PMF of beta-binomial with n=20, alpha=3, and beta=2 already known
pmf <- function(x) {
  value <- (x+2) * (x+1) * (21-x) / (2 * 23 * 22 * 21)
  return(value)
}
lines(range, pmf(range), type = "l", col = "red")   # Theoretical PMF
legend("topleft", c("Sample density", "Theoretical density"), lty = 1, col = c("black", "red"))
(theo_mean <- n*alpha / (alpha + beta))   # Theoretical mean
(theo_var <- n*alpha*beta * (alpha + beta + n) / ((alpha + beta)^2 * (alpha + beta + 1)))   # Theoretical variance


#### 2 ####

# (b)

# Calculate von Mises PDF for a given x, theta_1, and theta_2
von_Mises <- function(x, theta_1, theta_2=0) {
  # theta_2 defaults to 0 if no value provided
  bessel <- integrate(function(y) {exp(theta_1 * cos(y)) / (2*pi)}, lower = 0, upper = 2*pi)
  value <- exp(theta_1 * cos(x-theta_2)) / (2*pi*bessel$value)
  return(value)
}

# Find maximum value of von Mises with theta_2 = 0
fmax <- function(theta_1) {
  bessel <- integrate(function(y) {exp(theta_1 * cos(y)) / (2*pi)}, lower = 0, upper = 2*pi)
  value <- exp(theta_1) / (2*pi*bessel$value)
  return(value)
}

# Generate n draws from von Mises using accept-reject algorithm
rvon_Mises <- function(n, m, theta_1, theta_2) {
  # n = number of draws
  # m = maximum value the function can take, acts as upper bound for U(0,m)
  output <- c()
  i <- 0
  while(i < n) {   # Keep going until n successful draws
    y <- runif(1, min = 0, max = 2*pi)   # Draw y from unif(0,2*pi)
    u <- runif(1, min = 0, max = m)   # Draw u from unif(0,m)
    vm <- von_Mises(y, theta_1, 0)   # Calculate PDF at y
    if (u < vm) {   # If u < pdf at y, then store that draw in output
      i <- i+1
      output[i] <- y
    }   # Otherwise ignore this draw
  }
  output <- (output + theta_2) %% (2*pi)   # Shift output by theta_2 because location family
  return(output)
}

# (c)

# Make 10000 draws from accept-reject von Mises
theta_1 <- 2
theta_2 <- pi/4
set.seed(1)
samp2 <- rvon_Mises(n = iter, m = fmax(theta_1), theta_1 = theta_1, theta_2 = theta_2)
hist(samp2, probability = TRUE)   # Sample histogram
lines(density(samp2))   # Sample density
density(samp2)$x[which(density(samp2)$y == max(density(samp2)$y))]   # Sample mode, should = theta_2
(theta_2)

# Compare with theoretical
range2 <- seq(from = 0, to = 2*pi, by = 0.005)  
vM <- von_Mises(x = range2, theta_1 = theta_1, theta_2 = theta_2)
lines(range2, vM, type = "l", col = "red")   # Theoretical PDF
legend("topright", c("Sample density", "Theoretical density"), lty = 1, col = c("black", "red"))
range2[which(vM == max(vM))]   # Theoretical mode, should = theta_2


#### 3 ####

# (a)

# Generate m draws from Wishart_d(sigma, n) distribution using Johnson algorithm
rwish <- function(m, n, sigma) {
  # m = number of required draws
  # n = degrees of freedom
  # sigma = covariance matrix
  
  d <- nrow(sigma)
  
  values <- list()
  
  for (reps in 1:m) {

    # Calculate values for T matrix in a d^2 length vector
    Tvals <- c()
    for (i in 1:d) {
      for (j in 1:d) {
        index <- (j-1)*d + i
        if (i>j) {
          Tvals[index] <- rnorm(n = 1, mean = 0, sd = 1)
        } else if (i==j) {
          Tvals[index] <- sqrt(rchisq(n = 1, df = n-i+1))
        } else {
          Tvals[index] <- 0
        }
      }
    }
    
    T <- matrix(data = Tvals, nrow = d, ncol = d, byrow = FALSE)   # Convert vector into dxd matrix
    
    A <- T %*% t(T)   # TT' is Wishart_d(I_d, n)
    
    R <- chol(sigma)   # R is upper triangular Choleski decomposition
    L <- t(R)   # L is lower triangular, transpose of R
    
    values[[reps]] <- L %*% A %*% R   # LAL' is Wishart_d(sigma, n)
  
  }
  
  return(values)
}

# (b)

# Make 1000 draws from Wishart
iter3 <- 1000
df <- 8
sigma <- matrix(data = c(1, -0.5, 0.5, -0.5, 1, -0.5, 0.5, -0.5, 1), 
                nrow = 3, ncol = 3, byrow = FALSE)
set.seed(1)
samp3 <- rwish(m = iter3, n = df, sigma = sigma)

# (c)

# Theoretical mean = df * sigma:
(theo_mean_w <- df*sigma)

# Sample mean:
samp_mean_w <- matrix(data = 0, nrow = 3, ncol = 3)
for (k in 1:iter3) {
  mat <- samp3[[k]]
  samp_mean_w <- samp_mean_w + mat/iter3
}
samp_mean_w

# Theoretical variance Var_ij = df * (sigma_ij^2 + sigma_ii * sigma_jj):
theo_var_w <- matrix(data = 0, nrow = 3, ncol = 3)
for (i in 1:3) {
  for (j in 1:3) {
    theo_var_w[i,j] <- df * (sigma[i,j]^2 + sigma[i,i] * sigma[j,j])
  }
}
theo_var_w

# Sample variance:
samp_var_w <- matrix(data = 0, nrow = 3, ncol = 3)
for (k in 1:iter3) {
  mat <- samp3[[k]]
  samp_var_w <- samp_var_w + (mat - samp_mean_w)^2/iter3
}
samp_var_w

