source("msnet_simu_fun2.R")
library(Rmpi)
library(snow)

cl=makeCluster(30)
clusterEvalQ(cl,{library(Matrix);source("senet.R");source("utils.R")})
clusterExport(cl,c("ctr","c1","df")) #clusterCall(cl,print,a+b)
clusterApplyLB(cl, 1:nrow(df), simu,ctr,c1,df)
stopCluster(cl)

save(ctr,c1,df,file = "./results/df.RData")