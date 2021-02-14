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

sum(nc_sum_vec)

#B. Gauss-Quadrature
#Legendre polynomial order 7 can be found here https://en.wikipedia.org/wiki/Legendre_polynomials
#order 7 
#P(x)= 1/16*(4298x^7-693*x^5+315*x^3-35*x)
legendre_7_vec <- c(-35/16,0,315/16,0,-693/16,0,4298/16)
polyroot(legendre_7_vec)

#D.  Gauss - Legendre
#P'(X)=1/16*(30086x^6-3465*x^4+945*x^2-35)
#Weights are given by: 
#w[i] = 2/((1-x[i]^2)*(P'(x[i]))^2), i = 1,...,n


#Note for Gauss-Quadrature need to do a change of variables
