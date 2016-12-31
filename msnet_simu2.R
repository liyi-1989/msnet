source("msnet_simu_fun.R")


cl=makeCluster(30)
clusterEvalQ(cl,{library(Matrix);source("senet.R")})
clusterExport(cl,c("paras","betastar","ctr","P","GL","Lambda","Lambda2","Lambda3","df")) #clusterCall(cl,print,a+b)
clusterApplyLB(cl, 1:nrow(df), simu, T, paras,ctr,betastar,Lambda2,"fixed",T,df,"2")
stopCluster(cl)

