library(Matrix)
library(Rmpi)
library(snow)
#library(parallel) 
source("senet.R")

paras=list(px=50, py=10, n=200, lambda1=1000, lambda2=1000, lambdaly=0, noise.signal = 0.25, noise.response = 10)
px=paras$px
betastar=sapply(1:px,twobump,px)
#betastar=sapply(1:px,onebump,px)
#betastar=sapply(1:px,fourblock,px)
#betastar=sin(2*pi*(1:px)/px)

ctr=list(L = 1, use.gram = TRUE, maxiter = 5000, tol = 1e-5, init = NULL, sigma = 0.9)

Lpx=crossprod(fd(paras$px)); Lpy=crossprod(fd(paras$py))
Lx=kronecker(Lpx,diag(rep(1,paras$py))); Ly=kronecker(Lpy,diag(rep(1,paras$px)))
P=matrix(0,paras$px*paras$py,paras$px*paras$py)
for(i in 1:paras$px){ for(j in 1:paras$py){ P[(i-1)*paras$py+j,(j-1)*paras$px+i]=1 } } # Permutation matrix
paras$P=P

#Lx=NULL;Ly=NULL;P=NULL;

# GL=crossprod( fd(paras$px*paras$py) ) 
# Lambda = GL
# Lambda2 = GL+ t(P)%*%GL%*%P
# Lambda=Ly
# Lambda2=Ly+t(P)%*%Lx%*%P
# Lambda3=t(P)%*%Lx%*%P
# L = Ll2(Xt, Lambda, lambda2 = paras$lambda2) ### fixed step-size
# L2 = Ll2(Xt, Lambda2, lambda2 = paras$lambda2) ### fixed step-size

simu=function(jobid=1,fixdata=T,paras,ctr,betastar,Lx=NULL,Ly=NULL,stepsize = "fixed",par=F,df=NULL,type=NULL){
  write.table(NULL,paste0("./results/working_job_",jobid,"_type_",type,".txt"))
  cat("1. Set Parameters ...\n")
  px=paras$px; py=paras$py; n=paras$n; lambda1=paras$lambda1; lambda2=paras$lambda2; lambdaly=paras$lambdaly; P=paras$P;  
  if(par){
    lambda1=df[jobid,"lambda1"]
    lambda2=df[jobid,"lambda2"]
    lambdaly=df[jobid,"lambdaly"]
    n=df[jobid,"n"]
  }  
   
  B=matrix(0,px,py); Y=matrix(0,n,py); X=matrix(0,n,px); X=list(py);
  noise.signal = paras$noise.signal; noise.response = paras$noise.response
  
  n_test=n; Y_test=matrix(0,n_test,py); X_test=list(py); 
    
  dataseed=ifelse(fixdata==T,1,jobid)
   
  cat("2. Generate Data ...\n")  
  for(i in 1:py){
    #B[,i]=betastar[(1:px-i)%%px+1]
    #B[,i]=betastar*2*(0.5+abs(i-py/2)/py) 
    if(i<py/2){B[,i]=betastar*(0.5+abs(i-py/2)/py)}else{B[,i]=-betastar*(0.5+abs(i-py/2)/py)}  
    #B[round(py/2):py,]=0  
    data_wide=simblock1d2(n = n, p=px, noise.signal = noise.signal, noise.response = noise.response, B[,i],dataseed)
    Y[,i]=data_wide$y; X[[i]]=data_wide$Xt
    data_wide_test=simblock1d2(n = n, p=px, noise.signal = noise.signal, noise.response = noise.response, B[,i],dataseed+0.1)
    Y_test[,i]=data_wide_test$y; X_test[[i]]=data_wide_test$Xt
  }
  # Block to long
  data=list(Xt=as.matrix(bdiag(X)), y=c(Y), beta = c(B)) #Xt=bdiag(lapply(1:py,function(dumy,a){a},X)) # same X # npy*pxpy

  if(is.null(Ly)){
    Lpy= cor(Y)
    Ly=kronecker(Lpy,diag(rep(1,px)))  
  }
  if(is.null(Lx)){
    Lx= cor(data$Xt)
    Lambda=lambdaly*Ly+(1-lambdaly)*Lx  
  }else{
    Lambda=lambdaly*Ly+(1-lambdaly)*t(P)%*%Lx%*%P
  }    
  ctr$L=ifelse(stepsize=="fixed",Ll2(data$Xt, Lambda, lambda2=lambda2),NULL)
    
  cat("3. Fit the model ...\n")
  out = senet_fista(data$Xt, data$y, Lambda, lossfun=l2loss, gradfun=l2grad, lambda1=lambda1, lambda2=lambda2, stepsize=stepsize, control = ctr)
  Bhat=matrix(out$beta,byrow = FALSE,nrow = px, ncol = py)
  
  Y_test_hat=as.matrix(bdiag(X_test))%*%(out$beta)
  E_test=sqrt(mean((c(Y_test)-Y_test_hat)^2))
  #cat(head(data$Xt))
  #slist=list(fit=out,input=list(B=B,px=px,py=py,n=n,Lambda=Lambda,data=data_w),output=list(Bhat=Bhat,jobid=jobid))
      slist=list(fit=out,input=list(B=B,px=px,py=py,n=n,Lambda=Lambda,lambda1=lambda1,lambda2=lambda2,lambdaly=lambdaly),output=list(Bhat=Bhat,jobid=jobid,rmse_test=E_test,Y_test=Y_test,Y_test_hat=Y_test_hat))
  fname=paste0("./results/n",n,"px",px,"py",py,"job",jobid,"T",type)
  save(slist,file = paste0(fname,".RData"))
  write.table(NULL,paste0("./results/finishing_job_",jobid,"_type_",type,".txt"))
  return()
}

#l1s=c(0,1,10,100,1000,10000)
#l2s=c(0,1,10,100,1000,10000)
l1s=1000
l2s=1000
lys=c(0,0.25,0.5,0.75,1) #(0:10)/10 #
ns= (1:5)*100 #c(200,400,600) 

df=NULL
for(i1 in 1:length(l1s)){
  for(i2 in 1:length(l2s)){
      for(i3 in 1:length(lys)){
          for(i4 in 1:length(ns)){
               df=rbind(df,c(i1,i2,i3,i4,l1s[i1],l2s[i2],lys[i3],ns[i4]))
          }
      }
  }
}
colnames(df)=c("id1","id2","id3","id4","lambda1","lambda2","lambdaly","n")
