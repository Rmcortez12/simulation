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
