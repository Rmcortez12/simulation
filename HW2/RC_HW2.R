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
#create solution vector
y_gq <-  1:length(legendre_roots)
for(i in 1:length(legendre_roots)){
  if(i==1){
    y_gq[i] = 2
  } else if(i == length(legendre_roots)||i == 2){
    y_gq[i] = 0
  }else if(i%%2==0){
    y_gq[i] = 2/(i+1)
  }else {
    y_gq[i] = 0
  }
}

#create weight matrix
weight_matrix_gq <- function(x){
  len_x <- length(x)
  z_mat <- matrix(1,ncol = len_x, nrow = len_x)
  
  for(i in 1:len_x){
    z_mat[len_x,i] <- x[i]
  }
  for(i in 3:len_x-1){
    for(j in 1:len_x){
      z_mat[i,j] <- x[j]^i
    }
  }
  for(i in 1:len_x){
    z_mat[len_x,i] <- x[i]^(2*len_x-1)
  }
  return(z_mat)
}

z<-weight_matrix_gq(x)

#solve for weights
weights_gq <- solve(z,matrix(y_gq))
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
print(gq_sum(10,0,x,fx,weights_gq))


#D.  Gauss - Legendre
#P'(X)=1/16*(30086x^6-3465*x^4+945*x^2-35)
#Weights are given by: 
#w[i] = 2/((1-x[i]^2)*(P'(x[i]))^2), i = 1,...,n


#Note for Gauss-Quadrature need to do a change of variables
