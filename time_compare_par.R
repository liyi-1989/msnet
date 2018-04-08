source("time_compare.R")

cl=makeCluster(nrow(df))
clusterEvalQ(cl,{library(Matrix);source("senet.R")})
clusterExport(cl,c("paras","ctr","df")) 
clusterApplyLB(cl, 1:nrow(df), simu, T, paras,ctr,"fixed",T,df,"b",1)
stopCluster(cl)



#simu(jobid=1,fixdata=T,paras=paras,ctr=ctr,stepsize = "fixed",par=F,df=NULL,type=NULL,nrep=1)