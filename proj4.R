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

newt <- function(theta, func, grad, hess, ..., tol, fscale, maxit, max.half, eps){
  
  
}














