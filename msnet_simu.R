library(Matrix)
source("senet.R")

paras=list(
  px=20,
  py=20,
  n=500,
  lambda1=1,
  lambda2=100,
  noise.signal = 0.25,
  noise.response = 10
)
px=paras$px
betastar=sapply(1:px,twobump,px)
#betastar=sin(2*pi*(1:px)/px)

ctr=list(L = 1,
         use.gram = TRUE,
         maxiter = 10000, 
         tol = 1e-8,
         init = NULL,
         sigma = 0.9)

GL=crossprod( fd(paras$px*paras$py) ) 
P=matrix(0,paras$px*paras$py,paras$px*paras$py)
for(i in 1:paras$px){ for(j in 1:paras$py){ P[(i-1)*paras$py+j,(j-1)*paras$px+i]=1 } } # Permutation matrix
Lambda = GL
Lambda2 = GL+ t(P)%*%GL%*%P
# L = Ll2(Xt, Lambda, lambda2 = paras$lambda2) ### fixed step-size
# L2 = Ll2(Xt, Lambda2, lambda2 = paras$lambda2) ### fixed step-size

simu=function(paras,ctr,betastar,Lambda,stepsize = "fixed"){
  cat("1. Set Parameters ...\n")
  px=paras$px; py=paras$py; n=paras$n; B=matrix(0,px,py); Y=matrix(0,n,py); X=matrix(0,n,px); X=list(py); 
  lambda1=paras$lambda1; lambda2=paras$lambda2; noise.signal = paras$noise.signal; noise.response = paras$noise.response
   
  cat("2. Generate Data ...\n")
  for(i in 1:py){
    B[,i]=betastar[(1:px-i)%%px+1]
    data=simblock1d(n = n, p=px, noise.signal = noise.signal, noise.response = noise.response, B[,i])
    Y[,i]=data$y; X[[i]]=data$Xt
  }
  # Block to long
  data=list(Xt=as.matrix(bdiag(X)), y=c(Y), beta = c(B)) #Xt=bdiag(lapply(1:py,function(dumy,a){a},X)) # same X # npy*pxpy
  ifelse(stepsize=="fixed",ctr$L=Ll2(data$Xt, Lambda, lambda2=lambda2),ctr$L=NULL)
  cat("3. Fit the model ...\n")
  out = senet_fista(data$Xt, data$y, Lambda, lossfun=l2loss, gradfun=l2grad, lambda1=lambda1, lambda2=lambda2, stepsize=stepsize, control = ctr)
  Bhat=matrix(out$beta,byrow = FALSE,nrow = px, ncol = py)
  
  return(list(fit=out,B=B,Bhat=Bhat,px=px,py=py,n=n,Lambda=Lambda))
}


res1=simu(paras,ctr,betastar,Lambda)
res2=simu(paras,ctr,betastar,Lambda2)

image(res1$B)

norm(res1$B-res1$Bhat,"F")
image(res1$Bhat)

norm(res2$B-res2$Bhat,"F")
image(res2$Bhat)

