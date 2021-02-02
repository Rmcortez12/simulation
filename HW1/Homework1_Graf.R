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
(x <- (-b + sqrt(b^2 - 4*a*c)) / (2*a))

polyroot(c(c[1],b,a[1]))
polyroot(c(c[2],b,a[2]))
polyroot(c(c[3],b,a[3]))

(x2 <- (-b + sqrt((b+2*a) * (b-2*a))) / (2*a))
(x3 <- (-b + sqrt(b+2*a) * sqrt(b-2*a)) / (2*a))
(x4 <- (-b + exp((1/2)*(log(b+2*a) + log(b-2*a)))) / (2*a))
(x5 <- (-b/(2*a)) + sqrt((b^2)/(4*a^2) - c/a))
(x6 <- (-b/(2*a)) + sqrt(b/(2*a) + sqrt(c/a)) * sqrt(b/(2*a) - sqrt(c/a)))
(x7 <- (-b/sqrt(b-2*a) + sqrt(b+2*a)) * (sqrt(b-2*a) / (2*a)))
(x8 <- (-b/(2*a)) + exp((1/2)*(log(b+2*a) + log(b-2*a)) - log(2*a)))

