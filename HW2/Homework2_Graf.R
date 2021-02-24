# Ricardo Cortez & Ben Graf
# STA 6133 - Simulation and Statistical Computing
# Homework 2
# Due 26 Feb 2021

pacman::p_load(pacman, statmod, lme4)

#### 1 ####

# (a)

n <- 7
a <- 0
b <- 10
xa <- seq(from = a, to = b, by = (b-a)/(n-1))
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

integrate(prob1afunc, lower = a, upper = b)

# (d)

prob1dorig <- function(x) { return( (x+1) / x^4 ) }   # Define function to be integrated
prob1dfunc <- function(t) { return( exp(t)*(t+2) / (t+1)^4 ) }   # Define function for Gauss-Laguerre

lag10 <- gauss.quad(n = 10, kind = "laguerre", alpha = 0)   #Get weights and nodes
lag100 <- gauss.quad(n = 100, kind = "laguerre", alpha = 0)   #Get weights and nodes

(Q1d_n10 <- sum(lag10$weights * prob1dfunc(lag10$nodes)))
(Q1d_n100 <- sum(lag100$weights * prob1dfunc(lag100$nodes)))

integrate(prob1dorig, lower = 1, upper = Inf)


#### 2 ####

#Import and set up data
setwd("/Users/Ben/Library/Mobile Documents/com~apple~CloudDocs/Documents/UTSA Master's/Semester 4/STA 6133 Simulation & Statistical Computing/Homework/Homework 2")
toenail <- read.csv("./toenail.txt", header = FALSE, sep = " ")
names(toenail) <- c("patient", "response", "treatment", "time")

model <- glmer(response ~ treatment*time + (1 | patient), data = toenail, family = "binomial")
summary(model)
#plot(model)

preds <- predict(model, newdata = toenail)
binpreds <- preds>=0.5
table(toenail$response, binpreds)
mean(toenail$response == binpreds)


#### 3 ####

# (a)

my_rand <- function(m, F) {
  # The function arguments must be m = the number of required draws,
  # and F, a continuous and strictly increasing cdf.
  
  
}

# (b)

#Use your function to simulate one random sample of size 10000 from the Gamma(3, 1) distribution, and one random sample of size 10000 from the Gamma(0.3,1) distribution. Assess the adequacy of your output by overlaying in the same plot the histogram of your output and the density of the distribution.

my_rand(m = 10000, F = pgamma())


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
  f <- ifelse(x<2*a | x>2*b, 0, ifelse(x<a+b, (x-2*a)/(b-a)^2, (2*b-x)/(b-a)^2) )
  return(f)
}

range <- seq(from = 2*a, to = 2*b, by = 0.01)
lines(x = range, y = tri_pdf(range, a, b), col = 4)
legend("topright", legend = c("Sample density", "Actual PDF"), lty = 1, col = c(2,4), )



