##########################################################################
#  R implementation of the structured elastic net
#  'Feature Selection Guided by Structural Information'
#
#  Slawski, M., zu Castell, W., Tutz, G.
#  The Annals of Applied Statistics, 4(2), 1056-1080 
#  
#  Optimization is based on FISTA:
#  A. Beck and M. Teboulle 
#   'A Fast Iterative Shrinkage-Thresholding Algorithm for
#    Linear Inverse Problems' 
#  SIAM J. on Imaging Sciences, 2(1), 2009
#  
#  Any loss functions which is continuously
#  differentiable can be run with the code template provided
#  here. Squared loss and logistic loss are already implemented.
#
#  Arguments: 
#  X --- design matrix, including penalized and unpenalized (such
#        as the constant term) features
#  y --- response
#  Lambda ---  structure matrix for the quadratic regularizer
#  lossfun --- a continuously differentiable loss function
#              s. the examples below of how to implement
#              them
#  gradfun --- a function that evaluates the gradient of the
#              loss function
#  lambda1 --- regularization parameter for the ell_1-term
#  lambda2 --- regularization parameter for the quadratic term
#  stepsize --- either 'fixed' or 'backtracking'.
#               If 'stepsize == fixed' an upper bound on the
#               Lipschitz constant of the gradient has to
#               be provided. (s. 'control' below)
#  control: a list of OPTIONAL additional arguments
#  weights: feature-wise weights for the ell_1-term.
#           Weights equal to zero lead to removal of the
#           ell_1-penalty for the respective feature.
#           Can also be used for an adaptive version
#           of the structured elastic net.
#  L: Lipschitz constant of the gradient.
#  use.gram: whether the Gram matrix t(X) %*% X
#            should be computed stored.
#  maxiter: an upper bound on the number of iterations.  
#  tol: If successive iterates differ less than tol
#       (w.r.t. ell_2-norm), the algorithm stops.
#  init: starting value. Default is the zero vector.
#  sigma: parameter for back-tracking.
#
#   Value:
#
#   A list containing the following elements. 
#   beta: The estimate. 
#   iter: number of iterations run
#   obj: objective values of the iterates
#   objbetahat: the objective value of the estimate betahat.
#   kkt: KKT optimality
#   tol: difference of the last two successive iterates
#        (w.r.t. ell_2-norm)
#    L: Lipschitz constant of the gradient or an estimate thereof
#       (if stepsize was set to 'backtracking'). 
#  (C) Martin Slawski, Feb. 2012
#  Please report bugs to martin.dot.slawski (at) googlemail.com
##########################################################################
library(MASS)

senet_fista <- function(X, y, Lambda, lossfun, gradfun, lambda1, lambda2,
                        stepsize = c("fixed", "backtracking"),
                        control = list(weights = 1,
                                       L = NULL,
                                       use.gram = TRUE,
                                       maxiter = NULL, 
                                       tol = 1e-8,
                                       init = NULL,
                                       sigma = 0.9)){
#### input checking
  p <- ncol(X)
  n <- nrow(X)
  if(length(y) != nrow(X))
    stop("Dimensions of 'X' and 'y' do not match \n")
  if(!all(dim(Lambda) == p))
    stop("Dimensions of X and Lambda do not match \n")
  
  if(!all(c(lambda1, lambda2) >= 0))
    stop("Regularization parameters lambda1, lambda2 have
          to be non-negative \n")
  
  #### un-pack
  stepsize <- match.arg(stepsize)
  
  weights <- control$weights
  
  ###
  
  if(is.null(weights)){ 
  
    weights <- rep(1, p) 
    
  }
  
  else{
  
    if(length(weights) == 1)
       weights <- rep(weights, p)
    else
       if(!(length(weights) == p))
          stop("The length of 'weights' has to be equal to the number of columns of X \n")
          
     if(any(weights < 0))
      stop("The entries of weights have to be non-negative \n")     
  }
  
  ###
  L <- control$L
  if(is.null(L)){
    if(stepsize == "fixed"){ 
      cat("Lipschitz constant L unspecified; 'stepsize' set to 'backtracking' \n")
      stepsize <- "backtracking"
    }
  }
  use.gram <- control$use.gram
  if(is.null(use.gram)){
    use.gram <- FALSE
  }
  else
    XtX <- crossprod(X)
  maxiter <- control$maxiter
  if(is.null(maxiter))
    maxiter <- Inf
    
  tol <- control$tol
  if(is.null(tol))
    tol <- 1e-8
  init <- control$init
  if(!is.null(init)){
    if(!length(init) == p){
      warning('Dimension of initial solution does not match \n')
      
    }
    else
      beta <- init
  }
  else
    beta <- rep(0, p) ### initalize as 0-vector. 
    
  sigma <- control$sigma
  if(is.null(sigma)){
    sigma <- 0.9
    
  }  
  gamma <- beta  ### for the extrapolation steps.
  eta <- X %*% gamma
  Lambdatimesgamma  <- lambda2 * (Lambda %*% gamma)
  objs <- lossfun() +  (t(gamma) %*% Lambdatimesgamma) + lambda1 * sum(abs(gamma) * weights) ### objective
  eps <- .Machine$double.eps ### machine precision
  ###
  Xty <- t(X) %*% y
  
  ### KKT optimality
  kktvec <- gradfun() +  2 * Lambdatimesgamma
  A <- abs(gamma) > eps
  if(sum(A) > 0)
    kktA <- max(abs(kktvec[A] - lambda1 * weights[A]))
  else
    kktA <- 0
  if(sum(!A) > 0)
    kktAc <- max(pmax(abs(kktvec[!A]) - lambda1 * weights[!A], 0))
  else
    kktAc <- 0
  kktopt <- max(kktA, kktAc)
  ###
  betahat <- gamma ### the iterate to be returned.
  ###
  k <- 0 ### iteration counter
  stopcond <- FALSE
  t <- 1 ### initialization of stepsize (if back-tracking is used)
  
  ### helper function: soft thresholding operator
  soft <- function(x, alpha) pmax(abs(x) - alpha, 0) * sign(x)
  
  ### ** End of initialization steps ** ###
  ### Main FISTA algorithm
  while(!stopcond){
#    if(k>10){
#      if( (k%%floor(maxiter/10)==0)){
#        cat("Iteration ",k,"...\n")
#      }
#    }

    ### evaluate gradient
    grad <- gradfun() + 2 * lambda2 * Lambda %*% gamma
    ###
     
    
    if(stepsize == "fixed"){
      betanew <- soft(gamma - grad/L, lambda1 * weights/ L)
      eta <- X %*% betanew
      lossbetanew <- lossfun()
      ####
      Lambdatimesbetanew <- lambda2 * (Lambda %*% betanew)
      kktvec <- gradfun() +  2 * Lambdatimesbetanew
      
      A <- abs(betanew) > eps
       if(sum(A) > 0)
         kktA <- max(abs(kktvec[A] - lambda1 * weights[A]))
       else
         kktA <- 0
      if(sum(!A) > 0)
        kktAc <- max(pmax(abs(kktvec[!A]) - lambda1 * weights[!A], 0))
      else
        kktAc <- 0
      kktopt <- c(kktopt, max(kktA, kktAc))
      
      ###
      objnew <- lossbetanew + t(betanew) %*% Lambdatimesbetanew + lambda1 * sum(abs(betanew) * weights)
      if(objnew < min(objs))
        betahat <- betanew
      objs <- c(objs, objnew)
      
    }
    else{ ### back-tracking
      lossgamma  <-  lossfun()
      flag <- 1
      
        m <- 0
        while(flag  && t > eps){
            t <- t * sigma^m
            betanew <- soft(gamma - t * grad, lambda1 * weights * t)
            eta <- X %*% betanew
            lossbetanew  <- lossfun() 
            if(lossbetanew < lossgamma + sum(grad * (betanew - gamma)) + 1/(2* t) * sum((betanew - gamma)^2))
                flag <- 0
            else
                m <- m + 1
            
        }
        ####
        Lambdatimesbetanew <- lambda2 * (Lambda %*% betanew)
        kktvec <- gradfun() +  2 * Lambdatimesbetanew
      
        A <- abs(betanew) > eps
      if(sum(A) > 0)
         kktA <- max(abs(kktvec[A] - lambda1 * weights[A]))
       else
         kktA <- 0
      if(sum(!A) > 0)
        kktAc <- max(pmax(abs(kktvec[!A]) - lambda1 * weights[!A], 0))
      else
        kktAc <- 0
      kktopt <- c(kktopt, max(kktA, kktAc))
        
        ####
        objnew <- lossbetanew + t(betanew) %*% Lambdatimesbetanew + lambda1 * sum(abs(betanew) * weights)
         if(objnew < min(objs))
           betahat <- betanew
      objs <- c(objs, objnew)
           
    }
    gamma <- betanew + k/(k + 3) * (betanew - beta)
    delta <- sqrt( sum((beta - betanew)^2))
    if( delta < tol){
      beta <- betanew
      stopcond <- TRUE 
    }
    beta <- betanew
    k <- k + 1
    if(k > maxiter)
      stopcond <- TRUE
   eta <- X %*% gamma
      
  }  
    
  ### Return output
  out <- list()
  out$beta <- betahat
  out$iter <- k
  out$obj <- objs
  out$objbetahat <- min(objs)
  out$kkt <- kktopt
  out$tol  <- delta
  if(stepsize == "fixed")
    out$L <- L
  else
   out$L <- 1/t
  return(out)
  ###
} 
###################################################################  
### Example 1: squared loss
  
l2loss <- function() sum((get("y", parent.frame(1)) - get("eta", parent.frame(1)))^2)
l2grad <- function(){
  if(get("use.gram", parent.frame(1))){
    2 * (get("XtX", parent.frame(1)) %*% get("beta", parent.frame(1)) - get("Xty", parent.frame(1)) )
  }   
  else{
    2 * (t(get("X", parent.frame(1))) %*% get("eta", parent.frame(1)) - get("Xty", parent.frame(1)))
 }   
}
### Example 2: logistic loss
logisticloss <- function(){
  - sum(get("y", parent.frame(1)) *  get("eta", parent.frame(1)) - log(1 + exp(get("eta", parent.frame(1)))))
}
  
logisticgrad <- function(beta, use.gram){
  mu <- plogis( get("eta", parent.frame(1)) )
  t(get("X", parent.frame(1))) %*% mu - get("Xty", parent.frame(1))
}
### ... add a continuously differentiable loss function
##      with Lipschitz-continuous gradient of your
###     choice here ...
###################################################################
### helper functions
### computation of Lipschitz constant for squared loss
Ll2 <- function(X, Lambda, lambda2){
  #2 * (max( svd(crossprod(X) + lambda2 * Lambda)$d )^2) # this svd method seems slow
  2*base::norm(crossprod(X)+ lambda2 * Lambda,"2")^2
  
}
### computation of Lipschitz function for logistic loss 
Llogistic <- function(X, Lambda, lambda2){
    ( 0.25 * max( svd(crossprod(X) + lambda2 * Lambda)$d )^2)
  
  
}
###################################################################
###################################################################
  
### 1st example: two bumps

#### 
#simblock1d <- function(n, p=100, noise.signal = 0.25, noise.response = 30, beta,dataseed=1,...){
#  set.seed(dataseed)
#  Xt <- matrix(ncol = p, nrow = n)
#  for(i in 1:n){
#    bi <- runif(5, 0, 5)
#    mi <- runif(5, 0, 2*pi)
#    Xt[i,]  <- sapply(1:p, function(z) sum(bi * sin(z*pi*(5-bi)/50-mi))) + rnorm(p, sd=noise.signal)
#  }
#  y <- Xt %*% beta + rnorm(n, sd=noise.response)
#  return(list(Xt = Xt, y=y, beta = beta))
#}
#                      
#                      
#fd <- function(p){
#  D <- matrix(nrow = p-1, ncol = p, data = 0)
#  for(i in 1:(p-1)){
#    D[i,i] <- 1
#    D[i,i+1] <- -1
#  }
#  D
#} 
#
#twobump=function(t,px){
#  (floor(0.2*px) < t && t < floor(0.4*px))*((-1)*(floor(0.3*px) - t)^2+px)/(2*px) +
#  (floor(0.6*px) < t && t < floor(0.8*px))*((floor(0.7*px) - t)^2-px)/(2*px)
#}
#
#onebump=function(t,px){
#  (floor(0.2*px) < t && t < floor(0.4*px))*((-1)*(floor(0.3*px) - t)^2+px)/(2*px)
#}
## a=sapply(1:px, function(t) (20 < t && t < 40)*((-1)*(30 - t)^2+100)/200 +
##          (60 < t && t < 80)*((70 - t)^2-100)/200)
#
#fourblock=function(t,px){
#  (floor(0.2*px) < t && t <= floor(0.3*px))*0.5+
#  (floor(0.3*px) < t && t <= floor(0.4*px))*1+
#  (floor(0.4*px) < t && t <= floor(0.5*px))*0.5+
#  (floor(0.5*px) < t && t <= floor(0.6*px))*0.25
#}
#
#simblock1d2 <- function(n, p=100, noise.signal = 0.25, noise.response = 30, beta,dataseed=1,...){
#  set.seed(dataseed)
#  Xt <- matrix(ncol = p, nrow = n)
#  for(i in 1:n){
#    Mp=rep(0,p)  
#    Ip=diag(rep(1,p))
#    Xt[i,] = mvrnorm(1,Mp,Ip)
#  }
#  y <- Xt %*% beta + rnorm(n, sd=noise.response)
#  return(list(Xt = Xt, y=y, beta = beta))
#}







