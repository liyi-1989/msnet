library(Matrix)
library(Rmpi)
library(snow)
library(parallel) 
source("senet.R")

paras=list(
  px=100,
  py=20,
  n=200,
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

Lpx=crossprod(fd(paras$px))
Lpy=crossprod(fd(paras$py))
Lx=kronecker(Lpx,diag(rep(1,paras$py)))
Ly=kronecker(Lpy,diag(rep(1,paras$px)))
GL=crossprod( fd(paras$px*paras$py) ) 
P=matrix(0,paras$px*paras$py,paras$px*paras$py)
for(i in 1:paras$px){ for(j in 1:paras$py){ P[(i-1)*paras$py+j,(j-1)*paras$px+i]=1 } } # Permutation matrix
# Lambda = GL
# Lambda2 = GL+ t(P)%*%GL%*%P
Lambda=Ly
Lambda2=Ly+t(P)%*%Lx%*%P
Lambda3=t(P)%*%Lx%*%P
# L = Ll2(Xt, Lambda, lambda2 = paras$lambda2) ### fixed step-size
# L2 = Ll2(Xt, Lambda2, lambda2 = paras$lambda2) ### fixed step-size

simu=function(jobid=1,fixdata=T,paras,ctr,betastar,Lambda,stepsize = "fixed",par=F,df=NULL,type=NULL){
  write.table(NULL,paste0("./results/working_job_",jobid,".txt"))
  cat("1. Set Parameters ...\n")
  px=paras$px; py=paras$py; n=paras$n; B=matrix(0,px,py); Y=matrix(0,n,py); X=matrix(0,n,px); X=list(py); 
  lambda1=paras$lambda1; lambda2=paras$lambda2; noise.signal = paras$noise.signal; noise.response = paras$noise.response
  
  if(par){
    lambda1=df[jobid,3]
    lambda2=df[jobid,4]
  }
  
  dataseed=ifelse(fixdata==T,1,jobid)
   
  cat("2. Generate Data ...\n")
  for(i in 1:py){
    B[,i]=betastar[(1:px-i)%%px+1]
    data_w=simblock1d(n = n, p=px, noise.signal = noise.signal, noise.response = noise.response, B[,i],dataseed)
    Y[,i]=data_w$y; X[[i]]=data_w$Xt
  }
  # Block to long
  data=list(Xt=as.matrix(bdiag(X)), y=c(Y), beta = c(B)) #Xt=bdiag(lapply(1:py,function(dumy,a){a},X)) # same X # npy*pxpy
  ctr$L=ifelse(stepsize=="fixed",Ll2(data$Xt, Lambda, lambda2=lambda2),NULL)
  cat("3. Fit the model ...\n")
  out = senet_fista(data$Xt, data$y, Lambda, lossfun=l2loss, gradfun=l2grad, lambda1=lambda1, lambda2=lambda2, stepsize=stepsize, control = ctr)
  Bhat=matrix(out$beta,byrow = FALSE,nrow = px, ncol = py)
  
  cat(head(data$Xt))
  
  #slist=list(fit=out,input=list(B=B,px=px,py=py,n=n,Lambda=Lambda,data=data_w),output=list(Bhat=Bhat,jobid=jobid))
  slist=list(fit=out,input=list(B=B,px=px,py=py,n=n,Lambda=Lambda),output=list(Bhat=Bhat,jobid=jobid))
  fname=paste0("./results/n",n,"px",px,"py",py,"job",jobid,"T",type)
  save(slist,file = paste0(fname,".RData"))
  write.table(NULL,paste0("./results/finishing_job_",jobid,".txt"))
  return()
}

l1s=c(0,1,10,100,1000,10000)
l2s=c(0,1,10,100,1000,10000)

df=NULL
for(i in 1:length(l1s)){
  for(j in 1:length(l2s)){
    df=rbind(df,c(i,j,l1s[i],l2s[j]))
  }
}
colnames(df)=c("id1","id2","lambda1","lambda2")