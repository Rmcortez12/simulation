# Ben Graf & Ricardo Cortez
# STA 6113 - Simulation and Statistical Computing
# Homework 4
# Due 16 Apr 2021

pacman::p_load(pacman, nimble, metRology, lawstat)

#### 1 ####

# (a)

mu_0 <- 1   # The 4 distributions all have mean 1
set.seed(1)
M <- 100000
type1err <- data.frame("norm" = rep(0,4), "chisq" = rep(0,4), "unif" = rep(0,4), 
                       "exp" = rep(0,4), "alpha" = rep(0,4), "n" = rep(0,4))
config <- 0   # Counts different alpha and n configurations
for (alpha in c(0.1, 0.05)) {
  for (n in c(30, 300)) {
    config <- config + 1
    print(paste0("Configuration ",config))
    type1err[config,]$alpha <- alpha
    type1err[config,]$n <- n
    (crit_val <- qt(1-alpha/2, df = n-1))   # Critical value for t
    data1 <- matrix(rnorm(n*M, mean = 1, sd = sqrt(2)), nrow = M, ncol = n)
    data2 <- matrix(rchisq(n*M, df = 1), nrow = M, ncol = n)
    data3 <- matrix(runif(n*M, min = 0, max = 2), nrow = M, ncol = n)
    data4 <- matrix(rexp(n*M, rate = 1/1), nrow = M, ncol = n)
    t_test_stat1 <- (apply(data1, 1, mean) - mu_0) / (apply(data1, 1, sd) / sqrt(n))
    t_test_stat2 <- (apply(data2, 1, mean) - mu_0) / (apply(data2, 1, sd) / sqrt(n))
    t_test_stat3 <- (apply(data3, 1, mean) - mu_0) / (apply(data3, 1, sd) / sqrt(n))
    t_test_stat4 <- (apply(data4, 1, mean) - mu_0) / (apply(data4, 1, sd) / sqrt(n))
    type1err[config,1] <- mean(abs(t_test_stat1) > crit_val)
    type1err[config,2] <- mean(abs(t_test_stat2) > crit_val)
    type1err[config,3] <- mean(abs(t_test_stat3) > crit_val)
    type1err[config,4] <- mean(abs(t_test_stat4) > crit_val)
  }    
}

type1err


# (b)

alpha <- 0.05
n <- 300
set.seed(1)
M <- 10000
mu_range <- seq(from = 0, to = 4, by = 0.05)   # Range of "true" mu to estimate power for
crit_val <- qt(1-alpha/2, df = n-1)   # Critical value for t
power <- data.frame("norm" = rep(0,length(mu_range)), "chisq" = rep(0,length(mu_range)), 
                    "unif" = rep(0,length(mu_range)), "exp" = rep(0,length(mu_range)), "mu" = mu_range)

for (j in 1:length(mu_range)) {
  # The 4 functions have to be shifted to the current mean=mu being evaluated
  data1 <- matrix(rnorm(n*M, mean = mu_range[j], sd = sqrt(2)), nrow = M, ncol = n)
  data2 <- matrix(rchisq(n*M, df = mu_range[j]), nrow = M, ncol = n)
  data3 <- matrix(runif(n*M, min = mu_range[j]-1, max = mu_range[j]+1), nrow = M, ncol = n)
  data4 <- matrix(rexp(n*M, rate = 1/mu_range[j]), nrow = M, ncol = n)
  t_test_stat1 <- (apply(data1, 1, mean) - mu_0) / (apply(data1, 1, sd)/sqrt(n))
  t_test_stat2 <- (apply(data2, 1, mean) - mu_0) / (apply(data2, 1, sd)/sqrt(n))
  t_test_stat3 <- (apply(data3, 1, mean) - mu_0) / (apply(data3, 1, sd)/sqrt(n))
  t_test_stat4 <- (apply(data4, 1, mean) - mu_0) / (apply(data4, 1, sd)/sqrt(n))
  power[j,1] <- mean(abs(t_test_stat1) > crit_val)
  power[j,2] <- mean(abs(t_test_stat2) > crit_val)
  power[j,3] <- mean(abs(t_test_stat3) > crit_val)
  power[j,4] <- mean(abs(t_test_stat4) > crit_val)
  if (j %% 8 == 0) {
    print(j)
  }
}

plot(x = power$mu, y = power$norm, type = "l", lty = 1, col = 1, 
     xlab = expression(mu), ylab = "Power", 
     main = "Estimating Power of t-test for 4 families of distributions")
lines(x = power$mu, y = power$chisq, type = "l", lty = 1, col = 2)
lines(x = power$mu, y = power$unif, type = "l", lty = 1, col = 3)
lines(x = power$mu, y = power$exp, type = "l", lty = 1, col = 4)
legend("bottomright", c("Normal", "Chi-squared", "Uniform", "Exponential"), lty = 1, col = 1:4)


#### 2 ####

n <- 30
set.seed(1)
M <- 100000
theta <- 0   # MSEs should be constant for theta, so picking 0 for convenience


# (a)

# Normal
data_a <- matrix(rnorm(n*M, mean = theta, sd = sqrt(3)), nrow = M, ncol = n)
# Double exponential (Dexp)
data_b <- matrix(rdexp(n*M, location = theta, scale = sqrt(3/2)), nrow = M, ncol = n)
# Scaled T - for this function, sd is defined as "Scale factor for the shifted, scaled distribution"
data_c <- matrix(rt.scaled(n*M, df = 3, mean = theta, sd = 1), nrow = M, ncol = n)

# Function to estimate theta 3 ways and calculate MSE for each
MSEvec.func <- function(data, theta) {
  # data: Matrix of random data with each row containing a size n sample.
  #       Has M such rows and will estimate theta across these M Monte Carlo runs.
  # theta: The true value of the parameter being estimated.
  
  est_mean <- apply(data, 1, mean) 
  est_median <- apply(data, 1, median)
  est_trimmean <- apply(data, 1, mean, trim = 0.1)
  MSE_mean <- mean((est_mean - theta)^2)
  MSE_median <- mean((est_median - theta)^2)
  MSE_trimmean <- mean((est_trimmean - theta)^2)
  
  # Returns vector of the MSE of 3 estimates of theta:  mean, median, and trimmed mean (alpha=0.1)
  return(c(MSE_mean, MSE_median, MSE_trimmean))
}

MSEdf <- as.data.frame(rbind(MSEvec.func(data_a, theta),
                             MSEvec.func(data_b, theta),
                             MSEvec.func(data_c, theta)))
names(MSEdf) <- c("MSE_mean", "MSE_median", "MSE_trimmean")
rownames(MSEdf) <- c("Normal", "DExp", "Scaled T")
MSEdf


# (b)

set.seed(1)
epsilon_range <- c(0.1, 0.3)
MSEmix <- matrix(nrow = length(epsilon_range), ncol = 3)
for (k in 1:length(epsilon_range)) {
  mixing <- rbinom(M, size = 1, prob = epsilon_range[k])   # First draw from mixing distribution
  data_d <- matrix(nrow = M, ncol = n)
  for (i in 1:M) {
    if (mixing[i]) {   # If mixing drew a 1, draw from Scaled T
      data_d[i,] <- rt.scaled(n, df = 3, mean = theta, sd = 1)
    } else {   # If mixing drew a 0, draw from Normal
      data_d[i,] <- rnorm(n, mean = theta, sd = sqrt(3))
    }
  }
  MSEmix[k,] <- MSEvec.func(data_d, theta)   # Calculate 3 MSEs using function from (a)
}
MSEmix <- as.data.frame(MSEmix)
names(MSEmix) <- c("MSE_mean", "MSE_median", "MSE_trimmean")
rownames(MSEmix) <- c("Eps_0.1", "Eps_0.3")
MSEmix


#### 3 ####

k <- 3
n <- c(10, 10, 20)
alpha <- 0.05
M <- 10000


# (a)

groups <- c(rep(1,n[1]), rep(2,n[2]), rep(3,n[3]))
bart.test <- function(vec, groups, alpha) {
  return(bartlett.test(vec ~ groups)$p.value < alpha)   # If p-value < alpha, reject null hypothesis
}
lev.test <- function(vec, groups, alpha) {
  return(levene.test(vec, groups, location = "mean")$p.value < alpha)   # If p-value < alpha, reject null hypothesis
}

set.seed(1)

results <- matrix(nrow = 3, ncol = 4)

# Normal
var_0 <- c(1,1,1)
dat1 <- matrix(rnorm(n[1]*M, mean = 0, sd = sqrt(var_0[1])), nrow = M, ncol = n[1])
dat2 <- matrix(rnorm(n[2]*M, mean = 0, sd = sqrt(var_0[2])), nrow = M, ncol = n[2])
dat3 <- matrix(rnorm(n[3]*M, mean = 0, sd = sqrt(var_0[3])), nrow = M, ncol = n[3])
comb_dat <- cbind(dat1, dat2, dat3)
results[1,1] <- mean(apply(comb_dat, 1, bart.test, groups, alpha))   # Should approximate alpha = 0.05
results[1,3] <- mean(apply(comb_dat, 1, lev.test, groups, alpha))   # Should approximate alpha = 0.05

var_0 <- c(10,10,10)
dat1 <- matrix(rnorm(n[1]*M, mean = 0, sd = sqrt(var_0[1])), nrow = M, ncol = n[1])
dat2 <- matrix(rnorm(n[2]*M, mean = 0, sd = sqrt(var_0[2])), nrow = M, ncol = n[2])
dat3 <- matrix(rnorm(n[3]*M, mean = 0, sd = sqrt(var_0[3])), nrow = M, ncol = n[3])
comb_dat <- cbind(dat1, dat2, dat3)
results[1,2] <- mean(apply(comb_dat, 1, bart.test, groups, alpha))   # Should approximate alpha = 0.05
results[1,4] <- mean(apply(comb_dat, 1, lev.test, groups, alpha))   # Should approximate alpha = 0.05

# Scaled T
var_0 <- c(1,1,1)
dat1 <- matrix(rt.scaled(n[1]*M, df = 3, mean = 0, sd = sqrt(var_0[1])), nrow = M, ncol = n[1])
dat2 <- matrix(rt.scaled(n[2]*M, df = 3, mean = 0, sd = sqrt(var_0[2])), nrow = M, ncol = n[2])
dat3 <- matrix(rt.scaled(n[3]*M, df = 3, mean = 0, sd = sqrt(var_0[3])), nrow = M, ncol = n[3])
comb_dat <- cbind(dat1, dat2, dat3)
results[2,1] <- mean(apply(comb_dat, 1, bart.test, groups, alpha))   # Should approximate alpha = 0.05
results[2,3] <- mean(apply(comb_dat, 1, lev.test, groups, alpha))   # Should approximate alpha = 0.05

var_0 <- c(10,10,10)
dat1 <- matrix(rt.scaled(n[1]*M, df = 3, mean = 0, sd = sqrt(var_0[1])), nrow = M, ncol = n[1])
dat2 <- matrix(rt.scaled(n[2]*M, df = 3, mean = 0, sd = sqrt(var_0[2])), nrow = M, ncol = n[2])
dat3 <- matrix(rt.scaled(n[3]*M, df = 3, mean = 0, sd = sqrt(var_0[3])), nrow = M, ncol = n[3])
comb_dat <- cbind(dat1, dat2, dat3)
results[2,2] <- mean(apply(comb_dat, 1, bart.test, groups, alpha))   # Should approximate alpha = 0.05
results[2,4] <- mean(apply(comb_dat, 1, lev.test, groups, alpha))   # Should approximate alpha = 0.05

# Exponential
var_0 <- c(1,1,1)
dat1 <- matrix(rexp(n[1]*M, rate = 1/sqrt(var_0[1])), nrow = M, ncol = n[1])
dat2 <- matrix(rexp(n[2]*M, rate = 1/sqrt(var_0[2])), nrow = M, ncol = n[2])
dat3 <- matrix(rexp(n[3]*M, rate = 1/sqrt(var_0[3])), nrow = M, ncol = n[3])
comb_dat <- cbind(dat1, dat2, dat3)
results[3,1] <- mean(apply(comb_dat, 1, bart.test, groups, alpha))   # Should approximate alpha = 0.05
results[3,3] <- mean(apply(comb_dat, 1, lev.test, groups, alpha))   # Should approximate alpha = 0.05

var_0 <- c(10,10,10)
dat1 <- matrix(rexp(n[1]*M, rate = 1/sqrt(var_0[1])), nrow = M, ncol = n[1])
dat2 <- matrix(rexp(n[2]*M, rate = 1/sqrt(var_0[2])), nrow = M, ncol = n[2])
dat3 <- matrix(rexp(n[3]*M, rate = 1/sqrt(var_0[3])), nrow = M, ncol = n[3])
comb_dat <- cbind(dat1, dat2, dat3)
results[3,2] <- mean(apply(comb_dat, 1, bart.test, groups, alpha))   # Should approximate alpha = 0.05
results[3,4] <- mean(apply(comb_dat, 1, lev.test, groups, alpha))   # Should approximate alpha = 0.05

# Results
results <- as.data.frame(results)
names(results) <- c("Bartlett (1,1,1)", "Bartlett (10,10,10)", "Levene (1,1,1)", "Levene (10,10,10)")
rownames(results) <- c("Normal", "Scaled T", "Exponential")
results


# (b)

set.seed(1)
power3 <- matrix(nrow = 3, ncol = 4)

# Normal
var_0 <- c(1,1,3)
dat1 <- matrix(rnorm(n[1]*M, mean = 0, sd = sqrt(var_0[1])), nrow = M, ncol = n[1])
dat2 <- matrix(rnorm(n[2]*M, mean = 0, sd = sqrt(var_0[2])), nrow = M, ncol = n[2])
dat3 <- matrix(rnorm(n[3]*M, mean = 0, sd = sqrt(var_0[3])), nrow = M, ncol = n[3])
comb_dat <- cbind(dat1, dat2, dat3)
power3[1,1] <- mean(apply(comb_dat, 1, bart.test, groups, alpha))   # Gives power
power3[1,3] <- mean(apply(comb_dat, 1, lev.test, groups, alpha))   # Gives power

var_0 <- c(1,2,6)
dat1 <- matrix(rnorm(n[1]*M, mean = 0, sd = sqrt(var_0[1])), nrow = M, ncol = n[1])
dat2 <- matrix(rnorm(n[2]*M, mean = 0, sd = sqrt(var_0[2])), nrow = M, ncol = n[2])
dat3 <- matrix(rnorm(n[3]*M, mean = 0, sd = sqrt(var_0[3])), nrow = M, ncol = n[3])
comb_dat <- cbind(dat1, dat2, dat3)
power3[1,2] <- mean(apply(comb_dat, 1, bart.test, groups, alpha))   # Gives power
power3[1,4] <- mean(apply(comb_dat, 1, lev.test, groups, alpha))   # Gives power

# Scaled T
var_0 <- c(1,1,3)
dat1 <- matrix(rt.scaled(n[1]*M, df = 3, mean = 0, sd = sqrt(var_0[1])), nrow = M, ncol = n[1])
dat2 <- matrix(rt.scaled(n[2]*M, df = 3, mean = 0, sd = sqrt(var_0[2])), nrow = M, ncol = n[2])
dat3 <- matrix(rt.scaled(n[3]*M, df = 3, mean = 0, sd = sqrt(var_0[3])), nrow = M, ncol = n[3])
comb_dat <- cbind(dat1, dat2, dat3)
power3[2,1] <- mean(apply(comb_dat, 1, bart.test, groups, alpha))   # Gives power
power3[2,3] <- mean(apply(comb_dat, 1, lev.test, groups, alpha))   # Gives power

var_0 <- c(1,2,6)
dat1 <- matrix(rt.scaled(n[1]*M, df = 3, mean = 0, sd = sqrt(var_0[1])), nrow = M, ncol = n[1])
dat2 <- matrix(rt.scaled(n[2]*M, df = 3, mean = 0, sd = sqrt(var_0[2])), nrow = M, ncol = n[2])
dat3 <- matrix(rt.scaled(n[3]*M, df = 3, mean = 0, sd = sqrt(var_0[3])), nrow = M, ncol = n[3])
comb_dat <- cbind(dat1, dat2, dat3)
power3[2,2] <- mean(apply(comb_dat, 1, bart.test, groups, alpha))   # Gives power
power3[2,4] <- mean(apply(comb_dat, 1, lev.test, groups, alpha))   # Gives power

# Exponential
var_0 <- c(1,1,3)
dat1 <- matrix(rexp(n[1]*M, rate = 1/sqrt(var_0[1])), nrow = M, ncol = n[1])
dat2 <- matrix(rexp(n[2]*M, rate = 1/sqrt(var_0[2])), nrow = M, ncol = n[2])
dat3 <- matrix(rexp(n[3]*M, rate = 1/sqrt(var_0[3])), nrow = M, ncol = n[3])
comb_dat <- cbind(dat1, dat2, dat3)
power3[3,1] <- mean(apply(comb_dat, 1, bart.test, groups, alpha))   # Gives power
power3[3,3] <- mean(apply(comb_dat, 1, lev.test, groups, alpha))   # Gives power

var_0 <- c(1,2,6)
dat1 <- matrix(rexp(n[1]*M, rate = 1/sqrt(var_0[1])), nrow = M, ncol = n[1])
dat2 <- matrix(rexp(n[2]*M, rate = 1/sqrt(var_0[2])), nrow = M, ncol = n[2])
dat3 <- matrix(rexp(n[3]*M, rate = 1/sqrt(var_0[3])), nrow = M, ncol = n[3])
comb_dat <- cbind(dat1, dat2, dat3)
power3[3,2] <- mean(apply(comb_dat, 1, bart.test, groups, alpha))   # Gives power
power3[3,4] <- mean(apply(comb_dat, 1, lev.test, groups, alpha))   # Gives power

# Results
power3 <- as.data.frame(power3)
names(power3) <- c("Bartlett (1,1,3)", "Bartlett (1,2,6)", "Levene (1,1,3)", "Levene (1,2,6)")
rownames(power3) <- c("Normal", "Scaled T", "Exponential")
power3


#### 4 ####

fail_times <- c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)
xbar <- mean(fail_times)   # Sample mean
#xsd <- sd(fail_times)   # Sample standard deviation
n <- length(fail_times)   # Sample size
M <- 100000   # Simulation size


# (a)

set.seed(1)
pboot <- matrix(rexp(n*M, rate = 1/xbar), nrow = M, ncol = n)   # Parametric bootstrap sample
T.pboot <- 1/apply(pboot, 1, mean)   # M draws from bootstrap distribution of g(mu)=1/mu
hist(T.pboot)
(bias_est <- mean(T.pboot) - (1/xbar))
(sd_est <- sd(T.pboot))


# (b)

alpha <- 0.1   # 90% confidence intervals

# Exact
(exact_CI <- c(exp(-100 / (xbar * qgamma(alpha/2, shape = n, scale = 1/n))), exp(-100 / (xbar * qgamma(1-alpha/2, shape = n, scale = 1/n)))))

# Asymptotic
(asymp_CI <- c(exp(-100 / (qnorm(alpha/2) * xbar / sqrt(n) + xbar)), exp(-100 / (qnorm(1-alpha/2) * xbar / sqrt(n) + xbar))))


# Bootstrap
#hist(pboot)
gt_100 <- function(x) {
  # x is a vector
  return(mean(x>100))   # Returns percentage of elements of x that are >100
}
gt_fail <- gt_100(fail_times)   # Apply estimator to original data
Tb.pboot <- apply(pboot, 1, gt_100)
p_x_gt_100_est <- mean(Tb.pboot)   # Bootstrap mean of P(X>100), "center" of CI
p_x_gt_100_sd <- sd(Tb.pboot)   # Bootstrap standard deviation of P(X>100)
v <- unname(quantile(Tb.pboot, probs = c(1-alpha/2, alpha/2)))
c(2*gt_fail - v[1], 2*gt_fail - v[2])


# (c)

npboot <- matrix(sample(fail_times, size = n*M, replace = TRUE), nrow = M, ncol = n)   # Non-parametric bootstrap sample
Tc.npboot <- apply(npboot, 1, gt_100)
mean(Tc.npboot) 

# Bootstrap-t
V <- apply(npboot, 1, sd)/sqrt(n)
u <- unname(quantile((Tc.npboot - gt_fail)/V, probs = c(1-alpha/2, alpha/2)))
c(gt_fail - u[1]*sd(Tc.npboot)/sqrt(n), gt_fail - u[2]*sd(Tc.npboot)/sqrt(n))
## WHAT SHOULD GO WHERE sd(fail_times) or sd(Tc.npboot) IS CURRENTLY?

# Percentile
unname(quantile(Tb.pboot, probs = c(alpha/2, 1-alpha/2)))

# BCa
w <- qnorm(mean(Tc.npboot < gt_100(fail_times)))
k <- pnorm(2*w - qnorm(c(1-alpha/2, alpha/2)))
quantile(Tc.npboot, probs = k)
TJ <- c()
for(i in 1:n) {
  TJ[i] <- gt_100(fail_times[-i])
}
TJm <- mean(TJ)
a <- (1/6)*sum((TJm - TJ)^3) / sum((TJm - TJ)^2)^1.5
z <- qnorm(c(alpha/2, 1-alpha/2))
l <- pnorm(w + (w + z)/(1 - a*(w + z)))
quantile(Tc.npboot, probs = l)


