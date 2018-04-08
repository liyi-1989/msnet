library(Matrix)
library(Rmpi)
library(snow)
library(glmnet)
library(ggplot2)
#library(parallel) 
source("senet.R")
source("time_compare_fun.R")

paras=list(px=50, py=20, n=200, lambda1=100, lambda2=100, lambdaly=0.5, 
           noise.signal = 0.25, noise.response = 1)
ctr=list(L = 1, use.gram = TRUE, maxiter = 5000, tol = 1e-5, init = NULL, sigma = 0.9)

# try once
T1=simu(jobid=1,fixdata=T,paras=paras,ctr=ctr,stepsize = "fixed",par=F,df=NULL,type=NULL,nrep=1)

mean(T1)
sd(T1)

# try many 

pxs=(1:5)*10
pys=(1:2)*10
ns=c(200,400,600) 
nr=1:10


df=NULL
for(i1 in 1:length(ns)){
  for(i2 in 1:length(pys)){
    for(i3 in 1:length(pxs)){
      for(i4 in 1:length(nr)){
        
            df=rbind(df,c(i1,i2,i3,i4,ns[i1],pys[i2],pxs[i3],nr[i4]))

      }
    }
  }
}
colnames(df)=c("id1","id2","id3","id4","n","py","px","nrep")
save(df,file = "./results_time/results_df.RData")

T1list=list()
dfres=as.data.frame(matrix(0,nrow(df),2))
for(i in 1:nrow(df)){
  print(i)
  paras=list(px=df[i,"px"], py=df[i,"py"], n=df[i,"n"], lambda1=100, lambda2=100, lambdaly=0.5, 
             noise.signal = 0.25, noise.response = 1)
  T1=simu(jobid=1,fixdata=T,paras=paras,ctr=ctr,stepsize = "fixed",par=F,df=NULL,type=NULL,nrep=1)
  T1list[[i]]=T1
  dfres[i,]=c(mean(T1),sd(T1))
  save(T1list,df,dfres,file="results_time/T1list.RData")
}

df=cbind(df,dfres)
dfplot=subset(df,py==10,select = c("n","px","V1","V2"))
dfplot[,"n"]=as.factor(dfplot[,"n"])

# ggplot(dfplot, aes(x=px, y=V1, colour=n)) + 
#   geom_errorbar(aes(ymin=V1-V2, ymax=V1+V2),width=1) +
#   geom_line() +
#   geom_point()+
#   xlab("px") +
#   ylab("Time (sec)") +
#   ggtitle("Time Comparison")+
#   theme(plot.title = element_text(hjust = 0.5))




