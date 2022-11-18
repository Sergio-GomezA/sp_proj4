install.packages("numDeriv")
library('numDeriv')

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

theta0 <- c(0,0)
newt(theta0,func = rb,grad = gb,k=10,max.half = 5)

theta0 <- c(0,Inf)
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


########################################################################################################################################

#New testing functions
booth <- function(th){
  (th[1]+2*th[2]-7)^2 + (2*th[1]+th[2]-5)^2
}
boothgrad <- function(th){
  grad(booth, th)
}

multinv <- function(th){
  1/(th[1]*th[2])
}
multinvgrad <- function(th){
  grad(multinv, th)
}

saddle<- function(th){
  th[1]^2-th[2]^2
}
saddlegrad <- function(th){
  grad(saddle, th)
}

# Warning 1: Objective or derivatives not finite at theta0
theta0 = c(5,4)
newt(theta0, multinv, multinvgrad)
#WORKS for initial values non-finite with theta0=c(0,0)

# Warning 2: Step fails to reduce objective after trying max.half tries
# WORKS, tried using rb, gb with k=1e07+

# Warning 3: Maxit reached without convergence
newt(theta = c(1,1), booth, boothgrad, maxit = 100)
#WORKS, normal iterations required = 2 and fails at 1
newt(theta=c(1,3), booth, boothgrad, max.half = 1e+7)

optim(par=c(1,3), fn= booth,gr = boothgrad, method = "Nelder-Mead")
# Warning 4: Hessian not positive definite at convergence 
rb(theta0, k=100000)








