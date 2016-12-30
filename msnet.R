library(Matrix)
source("senet.R")

# Block Version
# Y: n*py
# X: n*px
# B: px*py
# n=dim(Y)[1]
# py=dim(Y)[2]
# px=dim(X)[2]

## 1. gnerate parameter
px=100
py=20
n=100
B=matrix(0,px,py)
Y=matrix(0,n,py)
X=matrix(0,n,px)
X=list(py)
P=matrix(0,px*py,px*py)

lambda1=10
lambda2=10
noise.signal = 0.25
noise.response = 10
ctr=list(L = L2,
     use.gram = TRUE,
     maxiter = 10000, 
     tol = 1e-8,
     init = NULL,
     sigma = 0.9)

GL=crossprod( fd(px*py) ) 
Lambda = GL
Lambda2 = GL+ t(P)%*%GL%*%P

L <- Ll2(Xt, Lambda, lambda2 = lambda2) ### fixed step-size
L2 <- Ll2(Xt, Lambda2, lambda2 = lambda2) ### fixed step-size

## 2. generate data: Y, X, with B

for(i in 1:px){
  for(j in 1:py){
    P[(i-1)*py+j,(j-1)*px+i]=1 # Permutation matrix 
  }
}

betastar <- sapply(1:px, function(t) (20 < t && t < 40)*((-1)*(30 - t)^2+100)/200 +
                     (60 < t && t < 80)*((70 - t)^2-100)/200)
for(i in 1:py){
  B[,i]=betastar[(1:px-i)%%px+1]
  data <- simblock1d(n = n, p=px, noise.signal = noise.signal, noise.response = noise.response, B[,i])
  Y[,i]=data$y
  X[[i]]=data$Xt
}

# Block to long
data=list(Xt=as.matrix(bdiag(X)), y=c(Y), beta = c(B)) #Xt=bdiag(lapply(1:py,function(dumy,a){a},X)) # same X # npy*pxpy

# 3. fit the model
#------------------------------
out <- senet_fista(data$Xt, data$y, Lambda, lossfun = l2loss, gradfun = l2grad, lambda1 = lambda1, lambda2 = lambda2,
                   stepsize = "fixed",
                   control = ctr)

Bhat=matrix(out$beta,byrow = FALSE,nrow = px, ncol = py)
#------------------------------

out2 <- senet_fista(data$Xt, data$y, Lambda2, lossfun = l2loss, gradfun = l2grad, lambda1 = lambda1, lambda2 = lambda2,
                   stepsize = "fixed",
                   control = ctr)

Bhat2=matrix(out2$beta,byrow = FALSE,nrow = px, ncol = py)

#------------------------------
image(B)

norm(B-Bhat,"F")
image(Bhat)

norm(B-Bhat2,"F")
image(Bhat2)



