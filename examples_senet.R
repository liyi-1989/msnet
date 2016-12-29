#########################################################################
#
# Examples for the R code implementing the structured elastic net
#
# Example 1: Two bumps
#
# Example 2: phonemes
#
# (C) Martin Slawski, Feb. 2012 
#
#########################################################################
#setwd("D:/works/senet")
source("senet.R")
### 1st example: two bumps
betastar <- sapply(1:100, function(t) (20 < t && t < 40)*((-1)*(30 - t)^2+100)/200 +
                                       (60 < t && t < 80)*((70 - t)^2-100)/200)

### 
set.seed(1809)
data <- simblock1d(n = 200, p=100, noise.signal = 0.25,
                   noise.response = 10, betastar)
X <- data$X
y <- data$y
Lambda <- crossprod( fd(100) )
### fixed step-size
L <- Ll2(X, Lambda, lambda2 = 10)
t0=proc.time()
out <- senet_fista(X, y, Lambda, lossfun = l2loss, gradfun = l2grad, lambda1 = 10^4, lambda2 = 10,
                        stepsize = "fixed",
                        control = list(L = L,
                                       use.gram = TRUE,
                                       maxiter = 10000, 
                                       tol = 1e-8,
                                       init = NULL,
                                       sigma = 0.9))
proc.time()-t0
t0=proc.time()
### back-tracking [recommended] 
out2 <- senet_fista(X, y, Lambda, lossfun = l2loss, gradfun = l2grad, lambda1 = 10^4, lambda2 = 10,
                        stepsize = "backtracking",
                        control = list(L = NULL,
                                       use.gram = TRUE,
                                       maxiter = 10000, 
                                       tol = 1e-8,
                                       init = NULL,
                                       sigma = 0.9))
proc.time()-t0
### 2nd example: phonemes
library(ElemStatLearn)
data(phoneme)
part <- phoneme[phoneme$g %in% c("ao", "aa"),]
X <- as.matrix(part[,1:256])
y  <- ifelse(part[,257] == "ao", 0, 1)
n <- length(y)
p <- ncol(X)
set.seed(1439)
trainind <- sample(1:n, floor(n*2/3))
Xtrain <- X[trainind,]
Xtest <-  X[-trainind,]
Ytrain <- y[trainind]
Ytest <- y[-trainind]
###
Lambda <- crossprod( fd(256) )
###
out2 <- senet_fista(Xtrain, Ytrain, Lambda, lossfun = logisticloss, gradfun = logisticgrad, lambda1 = 300, lambda2 = 25000,
                        stepsize = "backtracking",
                        control = list(L = NULL,
                                       use.gram = TRUE,
                                       maxiter = 10000, 
                                       tol = 1e-8,
                                       init = NULL,
                                       sigma = 0.9))
Yhat <-  (sign(Xtest %*% out2$beta)  + 1)/2
mean( abs(Yhat - Ytest) )
### add intercept (not penalized)
Lambdaprime <- rbind(0, cbind(0, Lambda))
myweights <- c(0, rep(1, p))
Xtrainprime <- cbind(1, Xtrain)
out2 <- senet_fista(Xtrainprime, Ytrain, Lambdaprime, lossfun = logisticloss, gradfun = logisticgrad, lambda1 = 300, lambda2 = 25000,
                        stepsize = "backtracking",
                        control = list(weights = myweights,
                                       L = NULL,
                                       use.gram = TRUE,
                                       maxiter = 10000, 
                                       tol = 1e-8,
                                       init = NULL,
                                       sigma = 0.9))
Yhat <-  (sign(cbind(1, Xtest) %*% out2$beta)  + 1)/2
mean( abs(Yhat - Ytest) )
#########################################################################
#########################################################################
 
 
 
 
