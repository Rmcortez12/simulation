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



