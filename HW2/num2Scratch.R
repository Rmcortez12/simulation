
# We are intested in approximating the log-likelihood using numerical integration methods

#Function to compute f(beta, y_i, sqrt(2)*sigma*t)
f<- function(t,eta,yi,trt,tx){
  #eta is a 3 length list with beta_0, beta_1,beta_2,beta_3,sigma
  #trt is treatment indicator [list] (294 elements)
  #t is time from the start of treatment of the jth visit for the ith patient is a jxi matrix (7 rows x 294 columns)
  a <- eta[1] + eta[2]*trt+eta[3]*trt*t+sqrt(2)*eta[4]*tx
  return(exp(sum(a*yi))/prod(1+a))
}

fv <- Vectorize(f, vectorize.args = "tx")

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
#patient ID, binary response, treatment code, time of eval
trt <- d[,3]
t <- matrix(nrow = 7, ncol = 294)
for(j in 1:7){
  for(i in 1:294){
    print(paste0("j: ",j))
    print(paste0("i: ",i))
    print(paste0("d: ",d[i,1]))
    if(d[i,1]==j){
      t[j,i]=d[i,4]
    }
  }
}

logliki_appx <- function(eta, yi){
  integral_appx <- sum(b$weights*fv(,eta,yi,trt,b$nodes))
  return(log(integral_appx))
}
