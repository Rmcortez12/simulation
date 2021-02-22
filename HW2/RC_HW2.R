#HW 2 Ricardo Cortez

### Number 1 ###
#approximate I using the Newton-Cotes quadrature rule of order n=7

#Find the quadrature points using 

fx <- function(x){
  #the function being integrated
  return(exp(-1/(1+x^2)))
}
true_integral <- integrate(fx,0,10)

#A. Newton-Cotes

quad_points_nc <- function(a,b,n){
  #a lower bound
  #b upper bound
  #n = order rule
  x<-1:n
  for(i in 1:n){
    x[i] <- a+(((i-1)*(b-a))/(n-1))
  }
  return(x)
}

weight_matrix_nc <- function(x){
  #x is the quadrature points
  len_x <- length(x)
  a <- x[1]
  b <- x[len_x]
  
  #create solution vector
  y <- 1:len_x
  for(i in 1:len_x){
    y[i] <- (b^i-a^i)/i
  }
  print(y)

  z_mat <- matrix(0, ncol=len_x+1, nrow = len_x)
  
  for(i in 1:length(z_mat[,1])){#rows
    print(i)
    for(j in 1:length(z_mat[1,])){#columns
      z_mat[i,j] <- x[j]^(i-1)
      z_mat[i,length(z_mat[1,])] <- y[i]
    }
  }
  return(z_mat)
}

x <- quad_points_nc(0,10,7)
z <- weight_matrix_nc(x)

weights <- solve(z[,1:7],z[,8])
fx_quads <- fx(x)
fx_quads
nc_sum_vec <- fx_quads*weights

nc_sum <-sum(nc_sum_vec)

#B. Gauss-Quadrature
#Legendre polynomial order 7 can be found here https://en.wikipedia.org/wiki/Legendre_polynomials
#order 7 
#P(x)= 1/16*(429x^7-693*x^5+315*x^3-35*x)
legendre_7_vec <- c(0,-35/16,0,315/16,0,-693/16,0,429/16)
legendre_roots <- Re(polyroot(legendre_7_vec))
x <- sort(legendre_roots)

#utilize the same matrix from part a just with legendre roots
xb_mat <- weight_matrix_nc(x)

#create solution vector
yb_vec <- 1:7
for(i in 1:7){
  yb_vec[i] <- (1^i-(-1)^i)/i
}

#solve for weights
wb_res <- solve(xb_mat[,1:7],yb_vec)



#weights_gq <- solve(z,matrix(y_gq))
#need to change intervals to fit [-1,1]
#f(x) --> f(5x+5)
gq_sum <- function(b,a,x,fx,w){
  #b upper bound
  #a lower bound
  #x regular quadratrue points
  #fx function
  #w weights
  # returns the estimate
  ba = (b-a)/2
  fx_transformed <- fx(((b-a)*x+a+b)/2)
  gq<- fx_transformed*w
  return(ba*sum(gq))
}


#SUmmary:
print("Integrate Function")
print(true_integral)
print("Newton-cotes:")
print(nc_sum)
print("Gauss-Quad")
print(gq_sum(10,0,x,fx,wb_res))


#D.  Gauss - Laguerre

#Define the equations
fx_d <- function(x){
  return((x+1)/x^4)
}
ft_d <- function(t){
  return(exp(t)*(t+2)/(t+1)^4)
}

#need points and weights for n=10, n=100
df10 <- gauss.quad(10, kind = "laguerre", alpha = 0)
df100 <- gauss.quad(100, kind = "laguerre", alpha = 0)

#now we sum
estimate10 <- sum(df10$weights*ft_d(df10$nodes))
estimate100 <- sum(df100$weights*ft_d(df100$nodes))

#Actual value
integrate(fx_d,1,Inf)
