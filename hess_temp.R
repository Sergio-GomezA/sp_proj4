

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

require(debug)
mtrace(newt)
newt(theta0,func = rb,grad = gb,hess = hb)
mtrace.off()


newt(theta0,func = rb,grad = gb,k=10,max.half = 5)



nll <- function(theta,t0,y) {
  ## -ve log likelihood for AIDS model y_i ~ Poi(alpha*exp(beta*t_i))
  ## theta = (alpha,beta)  
  mu <- theta[1] * exp(theta[2] * t0) ## mu = E(y)
  -sum(dpois(y,mu,log=TRUE)) ## the negative log likelihood
} ## nll

t80 <- 1:13 ## years since 1980
y <- c(12,14,33,50,67,74,123,141,165,204,253,246,240) ## AIDS cases

gll <- function(theta,t0,y) {
  ## grad of -ve log lik of Poisson AIDS early epidemic model
  alpha <- theta[1];beta <- theta[2] ## enhances readability
  ebt <- exp(beta*t0) ## avoid computing twice
  -c(sum(y)/alpha - sum(ebt),     ## -dl/dalpha
     sum(y*t0) - alpha*sum(t0*ebt)) ## -dl/dbeta
} ## gll

theta0 <- c(.5,0.5)
newt(theta0,func = nll,grad = gll,t0=t80,y=y)
