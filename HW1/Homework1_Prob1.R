#1 taylor expansion 

#Function for taylor expansion of sin(x)
taylor_expansion <- function(n,x){
  p <- (-1)^(n+1)*((x^(2*n-1))/factorial(2*n-1))
  return(p)
}

#Function for summation of expansion
loop_func <- function(n,x){
  sum = 0 
  for(i in 1:n){
    p = taylor_expansion(i,x)
    sum = sum+p
  }
  return(sum)
}

#Trying to consolidate data into data frame

sol <- data.frame("N"=c(10,20,100,10,20,100),"X"=c(pi,5*pi,pi,5*pi,pi,5*pi))

#x = pi
loop_func(10,pi)
loop_func(20,pi)
loop_func(100,pi)

#x = 5*pi
loop_func(10,5*pi)
loop_func(20,5*pi)
loop_func(100,5*pi)