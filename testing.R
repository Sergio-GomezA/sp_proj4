##################################################################################################################
################################################ STUFF TO TEST!!! ################################################
##################################################################################################################


## THE FOLLOWING SHOULD TEST IF A MATRIXIS A SQUARE MATRIX AND YIELD AN ERROR IF IT ISN'T
if(length(hess(theta)[1, ]) != length(hess(theta)[, 1])){
  stop("Hessian supplied is not a square matrix")
}

## THE FOLLOWING SHOULD TEST IF A HESSIAN IS SUPPLIED OR NOT 
## (comented out because the syntax yields an error w/o curly brackets)
#if(is.na(hess)==TRUE)
  
## THE FOLLOWING SHOULD TEST IF A HESSIAN IS OF FULL RANK, AND IT SHOULD PERTURB A MATRIX TO BE FULL
## RANK IF IT ISN'T

if(attr(chol(hessian, pivot=TRUE), "rank") != length(hessian[1, ])){
  
  ## Find a number which (when multiplied by the identity matrix)
  ## will make the eigen values of hess all greater than 0 by finding the min of the eigenvalues
  ## take the absolute value and add 1
  number <- abs(min(eigen(hessian)$values)) + 1
  
  ## construct an identity matrix of appropriate size
  iden <- diag(nrow=length(hessian[1,]))
  
  ## perturb hess to be positive definite by adding number*iden to it
  hessian <- hessian + number*iden
}

## THE FOLLOWING SHOULD GIVE THE INVERSE OF THE MATRIX HESSIAN
R <- chol(hessian)
inv_hess <- backsolve(R, forwardsolve(t(R), diag(nrow=length(hessian[1,]))))

## Should check that the delta is correct here
delta <- -inv_hess%*%grad

## Should check that theta + delta decreases func, if it does not halve it and check it again up to max.half times
## This shouldn't run at all if theta+delta decreases func
counter <- 0
while(func(theta+delta) >= func(theta)){
  delta <- delta/2
  counter <- counter+1
  if(counter == max.half){
    ## Not sure if I can break and then do stop but if I can this is how this should be because 
    ## I want the function to stop in that case
    break
    stop("max.half attempts done on current delta without decreasing the objective function")
  }
}














