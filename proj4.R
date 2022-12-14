## Em Belanger s2409525, Andreas Christoforou s2450189, Sergio Gomez s2441782
## github repo address: https://github.com/Sergio-GomezA/sp_proj4
## Em Belanger: Comments for newt function, detailed outline for newt function, 
## checks and warnings, perturb hessian, step halving when necessary, 
## create the iteration, convergence check, return
## Sergio Gomez: approx.Hess function and its testing. Testing and debugging 
## for newt function, testing extra "..." arguments. Using other functions for testing.
## Andreas Christoforou: Testing for errors with multiple new functions, 
## adding/correcting certain sections after doing the tests.


approx.Hess <- function(theta0,grad, eps = 1e-7,...){
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

################################################################################
## A function for minimizing functions using Newton's method.
##
## Inputs:
## theta: a vector of initial parameter values
## func: the objective function to minimize
## grad: the gradient function with the same arguments are func, but returns the 
##      gradient of the objective w.r.t.the elements of the parameter vector
## hess: (optional) the hessian matrix function, with the same arguments as func. 
##      If no hessian supplied the hessian will be estimated using finite 
##      differencing of the gradient vector (must be a square matrix)
## ... : any further arguments of func are passed using this
## tol: the convergence tolerance
## fscale: A rough estimate of the magnitude of func near the optimum used for 
##        convergence testing
## maxit: the maximum number of iterations to try before giving up
## max.half: the maximum number of halves to be taken from a step before 
##         concluding that the step has failed to improve the objective
## eps: the finite difference intervals to use when hess not supplied
## 
## Outputs:
## f: the value of the objective function at the minimum
## theta: a vector containing the values of the parameters at the minimum
## iter: the number of iterations taken to achieve the minimum
## g: the value of the gradient at the minimum
## Hi: the inverse of the hessian at the minimum
################################################################################



newt <- function(theta, func, grad, hess=NULL, ..., tol=1e-8, fscale=1, 
                 maxit=100, max.half=20, eps=1e-6){
  
  
  ## Check to see if the function or the gradient is finite at the initial theta,
  ## stop code if it is not finite
  
  if ((sum(is.finite(func(theta, ...))) != length(func(theta, ...))) 
      || (sum(is.finite(grad(theta, ...))) != length(grad(theta, ...)))){
    stop("The function or gradient is not finite at the initial theta.")
  }
  
  ## set iter = 1 to start (the first time running through the program is the 
  ## first iteration)
  iter <- 1
  
  ## Check if a hessian is supplied, and estimate a hessian if not
  if(is.null(hess)){
    hessian <- approx.Hess(theta, grad = grad, eps = eps,...)
  } 
  else{
    ## create the hessian matrix evaluated at theta as the starting hessian
    hessian <- hess(theta, ...)
    ## Check if hessian is square before we get started (and stop the program if
    ## it isn't with an error message)
    if(length(hess(theta)[1, ]) != length(hess(theta)[, 1])){
      stop("Hessian supplied is not a square matrix")
    }
  }
  
  ## Checks to see if the initial theta parameters are the minimum of the function,
  ## and if so returns the required values before the
  ## iteration process starts
  if(max(abs(grad(theta, ...))) < (tol*abs(func(theta, ...)+fscale))) {
    ## Evaluate the function at theta
    f <- func(theta, ...)
    ## Evaluate the gradient at theta
    g <- grad(theta, ...)
    
    #Checks to see whether the hessian is positive definite, assigns NA if not
    if(attr(chol(hessian, pivot=TRUE), "rank") != length(hessian[1, ])){
      warning("The hessian is not positive definite at the minimum")
      Hi <- NA
    }
    ## If the hessian is positive definite  calculate its inverse to be returned
    ## by the function
    else{
      Hi <- chol2inv(chol(hessian))
      ## Return f, theta, iter, g, and Hi for the initial theta, with iterations=0 
      ## as no iterations have occurred
      return (list("f" = f, "theta"=theta, "iter"= 0, "g"=g, "Hi"=Hi))
    }
  }
  
  ## iterate while iter <= maxit, break if max number of attempts reached
  while(iter <= maxit){
    
    ## Check that the hessian is positive definite by attempting a Cholesky Decomposition
    ## If it isn't positive definite perturb it to be so by multiplying by a 
    ## multiple of the identity matrix
    ## (sufficiently large multiple, check with Cholesky decomposition)
    
    ## using chol with pivot returns the rank of the cholesky decomposition 
    ## matrix as an attribute
    ## If the rank is not equal to the length of 1 row of the matrix the hessian
    ## is not of full rank
    ## Check this and perturb it to be pos def if it fails the test
    
    if(attr(chol(hessian, pivot=TRUE), "rank") != length(hessian[1, ])){
      
      ## Find a number which (when multiplied by the identity matrix)
      ## will make the eigen values of hess all greater than 0 by finding the 
      ## min of the eigenvalues
      ## take the absolute value and add 1
      number <- abs(min(eigen(hessian)$values)) + 1
      
      ## construct an identity matrix of appropriate size
      iden <- diag(nrow=length(hessian[1,]))
      
      ## perturb hess to be positive definite by adding number*iden to it
      hessian <- hessian + number*iden
    }
    
    ## Evaluate the expression: delta = negative inverse of the hessian multiplied 
    ## by the gradient
    ## We aren't using solve, so we will get the inverse of the hessian using cholesky 
    
    inv_hess <- chol2inv(chol(hessian))
    delta <- -inv_hess%*%grad(theta, ...)
    
    
    ## Check that theta + delta decreases func, if it does not halve it and 
    ## check it again up to max.half times
    counter <- 0
    fcurr <- func(theta, ...)
    fstep <- func(theta+delta, ...)
    if (!is.finite(fstep))
      stop("f is non finite at current step")
    while(is.finite(fstep) && fstep >= fcurr){
      delta <- delta/2
      counter <- counter+1
      if(counter == max.half){
        stop("max.half attempts done on current delta without decreasing the objective function")
      }
      fstep <- func(theta+delta, ...)
    }
    
    ## Update theta to be theta + delta
    theta <- theta + delta
    ## Update the hessian
    if(is.null(hess)){
      hessian <- approx.Hess(theta, grad = grad, eps = eps,...)
    }
    else{
      hessian <- hess(theta, ...)
    }
    
    ## Check if convergence reached by checking if all elements of the gradient 
    ## vector have absolute value
    ## less than tol*the absolute value of the objective function + fscale
    ## If convergence reached break and return the stuff
    
    if(max(abs(grad(theta, ...))) < (tol*abs(func(theta, ...)+fscale))){
      gg=grad(theta, ...)
      # print(paste("convergence reach at specified tolerance level.
      # iterations=",iter," gradient is ",gg[1],gg[2]))
      break
    }
    ## If convergence not reached increase iter by 1 and go through loop again
    else{
      iter <- iter + 1
    }
    
  }
  ## if convergence didn't happen and we iterated maxit times issue warning,
  ## but still return whatever we've
  ## managed to calculate
  if(iter>maxit){
    warning("Convergence not achieved after itertaing maximum number of times")
  }
  
  
  ## AFTER THE LOOP WHEN THE PROGRAM IS OVER ##
  ## check if the hessian is positive definite at the optimum and issue warning if it isn't
  ## If the matrix is not pos def at the minimum Hi gets "NA"
  if(attr(chol(hessian, pivot=TRUE), "rank") != length(hessian[1, ])){
    warning("The hessian is not positive definite at the minimum")
    Hi <- NA
  }
  ## If the hessian is positive definite at the minimum calculate its inverse 
  ## to be returned by the function
  else{
    Hi <- chol2inv(chol(hessian))
  }
  
  ## Evaluate the function at the minimum
  f <- func(theta, ...)
  
  ## Evaluate the gradient at the optimum
  g <- grad(theta, ...)
  
  ## Return f, theta, iter, g, and Hi
  return (list("f" = f, "theta"=theta, "iter"=iter, "g"=g, "Hi"=Hi))
  
}

