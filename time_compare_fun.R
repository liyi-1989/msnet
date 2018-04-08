simu=function(jobid=1,fixdata=T,paras,ctr,stepsize = "fixed",par=F,df=NULL,type=NULL,nrep=10){
  write.table(NULL,paste0("./results_time/working_job_",jobid,"_type_",type,".txt"))
  cat("1. Set Parameters ...\n")
  px=paras$px; py=paras$py; n=paras$n; lambda1=paras$lambda1; lambda2=paras$lambda2; lambdaly=paras$lambdaly; P=paras$P; 
  px=df[jobid,"px"]; py=df[jobid,"py"]; n=df[jobid,"n"]
  # if(par){
  #   lambda1=df[jobid,"lambda1"]
  #   lambda2=df[jobid,"lambda2"]
  #   lambdaly=df[jobid,"lambdaly"]
  #   n=df[jobid,"n"]
  #   px=py=df[jobid,"p"]
  #   dataseed=df[jobid,"dataseeds"]
  # }  
  dataseed=1
  B=matrix(0,px,py); Y=matrix(0,n,py); X=matrix(0,n,px); X=list(py);
  noise.signal = paras$noise.signal; noise.response = paras$noise.response
  
  n_test=n; Y_test=matrix(0,n_test,py); X_test=list(py); 
  
  #dataseed=ifelse(fixdata==T,1,jobid)
  
  
  Lpx=crossprod(fd(px)); Lpy=crossprod(fd(py))
  Lx=kronecker(Lpx,diag(rep(1,py))); Ly=kronecker(Lpy,diag(rep(1,px)))
  P=matrix(0,px*py,px*py)
  for(i in 1:px){ for(j in 1:py){ P[(i-1)*py+j,(j-1)*px+i]=1 } } # Permutation matrix
  
  cat("2. Generate Data ...\n")  
  
  betastar=sapply(1:px,onebump,px)
  
  for(i in 1:py){
    B[,i]=betastar[(1:px-i)%%px+1]
    #B[,i]=betastar*2*(0.5+abs(i-py/2)/py) 
    #if(i<py/2){B[,i]=betastar*(0.5+abs(i-py/2)/py)}else{B[,i]=-betastar*(0.5+abs(i-py/2)/py)}  
    #B[round(py/2):py,]=0  
    data_wide=simblock1d2(n = n, p=px, noise.signal = noise.signal, noise.response = noise.response, B[,i],dataseed)
    Y[,i]=data_wide$y; X[[i]]=data_wide$Xt
    # data_wide_test=simblock1d2(n = n, p=px, noise.signal = noise.signal, noise.response = noise.response, B[,i],dataseed+0.1)
    # Y_test[,i]=data_wide_test$y; X_test[[i]]=data_wide_test$Xt
  }
  # Block to long
  data=list(Xt=as.matrix(bdiag(X)), y=c(Y), beta = c(B)) #Xt=bdiag(lapply(1:py,function(dumy,a){a},X)) # same X # npy*pxpy
  
  #=================
  # Bhatlasso=matrix(NA,px,py)
  # for(j in 1:py){
  #   cvfit=cv.glmnet(X[[j]], Y[,j],nfolds = 5,alpha=1)
  #   #print(cvfit$lambda.1se)
  #   Bhatlasso[,j]=coefficients(cvfit)[-1]
  # }
  #=================
  lambda1=lambda2=1#cvfit$lambda.min
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
  
  T1=rep(0,nrep)
  for(i in 1:nrep){
    t0=proc.time()
    out = senet_fista(data$Xt, data$y, Lambda, lossfun=l2loss, gradfun=l2grad, lambda1=lambda1, lambda2=lambda2, stepsize=stepsize, control = ctr)
    tmp=proc.time()-t0
    T1[i]=tmp[3]
  }
  
  
  Bhat=matrix(out$beta,byrow = FALSE,nrow = px, ncol = py)
  
  # Y_test_hat=as.matrix(bdiag(X_test))%*%(out$beta)
  # E_test=sqrt(mean((c(Y_test)-Y_test_hat)^2))
  slist=list(time=T1,fit=out,input=list(B=B,px=px,py=py,n=n,Lambda=Lambda,lambda1=lambda1,lambda2=lambda2,lambdaly=lambdaly),output=list(Bhat=Bhat,jobid=jobid))
  fname=paste0("./results_time/n",n,"px",px,"py",py,"seed",dataseed,"job",jobid)
  save(slist,file = paste0(fname,".RData"))
  write.table(NULL,paste0("./results_time/finishing_job_",jobid,".txt"))
  
  print(paste("iter",out$iter))
  return(T1)
}

