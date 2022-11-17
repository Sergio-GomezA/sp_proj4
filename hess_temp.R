

approx.Hess <- function(theta0,grad,...){
  #'
  #'Calculates a numerical approximation of Hessian at value "theta"
  #'Using the gradient "grad"
  #'
  #'Inputs:
  #'theta0 : the point "vector" where the Hessian is going to be approximated 
  #'grad : a function that can evaluate the gradient of a function
  #'... : any other parameters that grad function might need
  #'
  #'Outputs:
  #'Hessy: Numerical approximation of the Hessian
  #'  
 
  eps <- 1e-7  ## finite difference interval 
  sizy <- length(theta0) # getting the dimensions
  gll0 <- grad(theta0,...) ## grad evaluation at theta0
  Hessy <- matrix(0,sizy,sizy) # initializing hessian
  for (i in 1:sizy) { ## loop over parameters
    th1 <- theta0; th1[i] <- th1[i] + eps ## increase theta0[i] by eps
    gll1 <- grad(th1,...) ## compute resulting nll
    Hessy[i,] <- (gll1 - gll0)/eps ## approximate second derivatives
  }
  
  return (Hessy)
}

rb <- function(th,k=2) {
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
}
gb <- function(th,k=2) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}
hb <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}


# testing

thets <- c(1,1)
rb(thets)
gb(thets)
hb(thets)


require(debug)
mtrace(approx.Hess)
mtrace.off()
approx.Hess(thets,gb)


#### missing
#### try in cholesky
#### non finite f evaluation !is.finite()



rb <- function(th,k=2) {
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
}
gb <- function(th,k=2) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}
hb <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}

theta0 <- c(0,0)

mtrace(newt)
newt(theta0,func = rb,grad = gb,hess = hb)
