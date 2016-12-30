library(Matrix)
source("senet.R")

paras=list(
  px=10,
  py=10,
  n=50,
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

simu=function(jobid=1,fixdata=T,paras,ctr,betastar,Lambda,stepsize = "fixed",par=F,df=NULL,type=NULL){
  
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
  
  slist=list(fit=out,input=list(B=B,px=px,py=py,n=n,Lambda=Lambda,data=data_w),output=list(Bhat=Bhat,jobid=jobid))
  fname=paste0("./results/n",n,"px",px,"py",py,"job",jobid,"T",type)
  save(slist,file = paste0(fname,".RData"))
  return(slist)
}



# Test serial
t0=proc.time()
for(i in 1:2){
  for(j in 1:2){
    cat("============== lambda 1:", i, "lambda 2:", j, " ==============\n")
    paras$lambda1=10^(i-1)
    paras$lambda2=10^(j-1)
    #res=simu(jobid=paste0(i,"_",j),F,paras,ctr,betastar,Lambda)
    res=simu(jobid=10*i+j,F,paras,ctr,betastar,Lambda)
  }
}
t1=proc.time()-t0



# Test parallel
library(Rmpi)
library(snow)
library(parallel) 

t0=proc.time()

l1s=c(1,10)
l2s=c(1,10)

df=NULL
for(i in 1:length(l1s)){
  for(j in 1:length(l2s)){
    df=rbind(df,c(i,j,l1s[i],l2s[j]))
  }
}
colnames(df)=c("id1","id2","lambda1","lambda2")

cl=makeCluster(3)
clusterEvalQ(cl,{library(Matrix);source("senet.R")})
clusterExport(cl,c("paras","betastar","ctr","P","GL","Lambda","Lambda2","df")) #clusterCall(cl,print,a+b)
clusterApplyLB(cl, 1:nrow(df), simu, F, paras,ctr,betastar,Lambda,"fixed",T,df,"1")
clusterApplyLB(cl, 1:nrow(df), simu, F, paras,ctr,betastar,Lambda2,"fixed",T,df,"2")
stopCluster(cl)

t2=proc.time()-t0


# Test code

# res1=simu(paras,ctr,betastar,Lambda)
# res2=simu(paras,ctr,betastar,Lambda2,jobid = 2)
# 
# image(res1$input$B)
# 
# cat(norm(res1$input$B-res1$output$Bhat,"F"),"\n")
# image(res1$output$Bhat)
# 
# cat(norm(res2$input$B-res2$output$Bhat,"F"),"\n")
# image(res2$output$Bhat)

# read in result

for(jobid in 1:nrow(df)){
  fname=paste0("./results/n",paras$n,"px",paras$px,"py",paras$py,"job",jobid,"T","2")
  load(paste0(fname,".RData"))
  cat("Job",jobid,":",norm(slist$input$B-slist$output$Bhat,"F"),"\n")
}

load("./results/n50px10py10job4T1.RData")
head(slist$input$data$Xt)