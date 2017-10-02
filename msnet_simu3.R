source("msnet_simu_fun.R")


cl=makeCluster(100)
clusterEvalQ(cl,{library(Matrix);source("senet.R")})
clusterExport(cl,c("paras","betastar","ctr","P","GL","Lambda","Lambda2","Lambda3","df")) #clusterCall(cl,print,a+b)
clusterApplyLB(cl, 1:nrow(df), simu, T, paras,ctr,betastar,Lambda3,"fixed",T,df,"3")
stopCluster(cl)

