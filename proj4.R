## NAMES
## github repo address: https://github.com/Sergio-GomezA/sp_proj4 (I THINK IT'S THIS BUT SOMEONE CONFIRM?!)
## CONTRIBUTIONS!


##################################################################################################################
## A function for minimizing functions using Newton's method.
##
## Inputs:
## theta: a vector of initial parameter values
## func: the objective function to minimize
## grad: the gradient function with the same arguments are func, but returns the gradient of the objective w.r.t.
## the elements of the parameter vector
## hess: (optional) the hessian matrix function, with the same arguments as func. If no hessian supplied the
## hessian will be estimated using finite differencing of the gradient vector
## ... : any further arguments of func are passed using this
## tol: the convergence tolerance
## fscale: A rough estimate of the magnitude of func near the optimum used for convergence testing
## maxit: the maximum number of iterations to try before giving up
## max.half: the maximum number of halves to be taken from a step before concluding that the step has failed 
## to improve the objective
## eps: the finite difference intervals to use when hess not supplied
## 
## Outputs:
## f: the value of the objective function at the minimum
## theta: a vector containing the values of the parameters at the minimum
## iter: the number of iterations taken to achieve the minimum
## Hi: the inverse of the hessian at the minimum
##################################################################################################################
## NOT FOR US TO DELETE LATER!!!
## Issue warning using "stop" or "warning" if objective or derivatives are not finite at the initial theta
## If step fails to reduce the function even after trying max.half step halvings
## If maxit is reached w/o convergence
## If the hessian is not positive definite at the convergence
##################################################################################################################

newt <- function(theta, func, grad, hess, ..., tol, fscale, maxit, max.half, eps){
  
  ## set iter = 0 to start
  iter <- 0
  
  ## Check if a hessian is supplied, and estimate a hessian if not
  if(is.na(hess)==TRUE){
    
  }
  ## In the final version this will be where iteration starts
  
  ## Check that the hessian is positive definite by attempting a Cholesky Decomposition
  ## If it isn't positive definite perturb it to be so by multiplying by a multiple of the identity matrix
  ## (sufficiently large multiple, check with Cholesky decomposition)
  
  ## using chol with pivot returns the rank of the cholesky decomposition matrix as an attribute
  ## If the rank is not equal to the length of 1 row of the matrix the hessian is not of full rank
  ## Check this and perturb it to be pos def if it fails the test
  
  if(attr(chol(rank, pivot=TRUE), "rank") != length(hess[1, ])){
    
    ## Find a number which (when multipled by the identity matrix)
    ## will make the eigen values of hess all greater than 0 by finding the min of the eigenvalues
    ## take the absolute value and add 1
    number <- abs(min(eigen(hess)$values, "values")) + 1
    
    ## construct an identity matrix of appropriate size
    iden <- diag(nrow=length(hess[1, ]))
    
    ## perturb hess to be positive definite by adding number*iden to it
    hess <- hess + number*iden
  }
  
  ## Evaluate the expression: delta = negative inverse of the hessian multiplied by the gradient
  ## We aren't using solve, so we will get the inverse of the hessian using cholesky and backsolve/forwardsolve
  R <- chol(test)
  inv_hess <- backsolve(R, forwardsolve(t(R), diag(nrow=length(hess[1,]))))
  delta <- -inv%*%grad
  
  
  ## Check that theta + delta decreases func, if it does not halve it and check it again up to max.half times
  
  ## Check if convergence reached by checking if all elements of the gradient vector have absolute value
  ## less than tol*the absolute value of the objective function + fscale
  
  ## If convergence reached break and return the stuff
  ## If convergence not reached increase iter by 1 and go through loop again
  
}













