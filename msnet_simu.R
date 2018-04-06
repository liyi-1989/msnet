source("msnet_simu_fun.R")
# Test parallel
#t0=proc.time()
cl=makeCluster(nrow(df))
clusterEvalQ(cl,{library(glmnet);library(Matrix);source("senet.R")})
clusterExport(cl,c("paras","ctr","df")) #clusterCall(cl,print,a+b)
clusterApplyLB(cl, 1:nrow(df), simu, T, paras,ctr,"fixed",T,df,"b")
#clusterApplyLB(cl, 1:nrow(df), simu, T, paras,ctr,betastar,Lx,Ly,"fixed",T,df,"b")
stopCluster(cl)
#t2=proc.time()-t0


# Test code
# res1=simu(jobid=1,fixdata=T,paras,ctr,stepsize = "fixed",par=T,df=df,type=NULL)
# #res2=simu(jobid=1,fixdata=T,paras,ctr,betastar,Lambda2)
# 
# par(mfrow=c(1,2))
# image(res1$input$B)
# image(res1$output$Bhat)
# plot(res1$input$B[1,],col="blue")
# plot(res1$output$Bhat[1,],col="red")
# printM(res1$input$B)
# printM(res1$output$Bhat)

# cat(norm(res1$input$B-res1$output$Bhat,"F"),"\n")
# image(res1$output$Bhat)
# 
# cat(norm(res2$input$B-res2$output$Bhat,"F"),"\n")
# image(res2$output$Bhat)

# read in result

#for(jobid in 1:nrow(df)){
#  fname=paste0("./results/n",paras$n,"px",paras$px,"py",paras$py,"job",jobid,"T","1")
#  load(paste0(fname,".RData"))
#  cat("Job",jobid,":",norm(slist$input$B-slist$output$Bhat,"F"),"\n")
#}
#
#
#for(jobid in 1:nrow(df)){
#  fname=paste0("./results/n",paras$n,"px",paras$px,"py",paras$py,"job",jobid,"T","2")
#  load(paste0(fname,".RData"))
#  cat("Job",jobid,":",norm(slist$input$B-slist$output$Bhat,"F"),"\n")
#}
#
#for(jobid in 1:nrow(df)){
#  fname=paste0("./results/n",paras$n,"px",paras$px,"py",paras$py,"job",jobid,"T","3")
#  load(paste0(fname,".RData"))
#  cat("Job",jobid,":",norm(slist$input$B-slist$output$Bhat,"F"),"\n")
#}

#load("./results/n50px10py10job4T1.RData")
#head(slist$input$data$Xt)
#load("./results/n50px10py10job4T2.RData")
#head(slist$input$data$Xt)
#load("./results/n50px10py10job3T2.RData")
#head(slist$input$data$Xt)
#load("./results/n50px10py10job3T1.RData")
#head(slist$input$data$Xt)

# Test serial
#t0=proc.time()
#for(i in 1:2){
#  for(j in 1:2){
#    cat("============== lambda 1:", i, "lambda 2:", j, " ==============\n")
#    paras$lambda1=10^(i-1)
#    paras$lambda2=10^(j-1)
#    #res=simu(jobid=paste0(i,"_",j),F,paras,ctr,betastar,Lambda)
#    res=simu(jobid=10*i+j,F,paras,ctr,betastar,Lambda)
#  }
#}
#t1=proc.time()-t0