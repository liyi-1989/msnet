library(Matrix)
library(glmnet)
library(ggplot2)
source("senet.R")

#load("../clinet/data/air_mon_mean_mon_mean_not_removed_sub.RData")
load("../clinet/data/air_mon_mean_mon_mean_removed_sub.RData")
load("../obs_data_mon_avg/monavg_tmmx_1979_2008.RData")
load("../obs_data_mon_avg/monavg_tmmx_2009_2016.RData")

plon=-71.0589 ;plat=42.3601; Plon=-70; Plat=40 ;CN="Boston"
plon=-71.4676 ;plat=42.7654; Plon=-70; Plat=45 ;CN="Nashua"
plon=-71.4128 ;plat=41.8240; Plon=-70; Plat=40 ;CN="Providence"
plon=-74.0060 ;plat=40.7128; Plon=-75; Plat=40 ;CN="New York" 
plon=-77.0369 ;plat=38.9072; Plon=-75; Plat=40 ;CN="Washington, D.C."
plon=-72.6851 ;plat=41.7637; Plon=-75; Plat=40 ;CN="Hartford"
plon=-118.2437 ;plat=34.0522; Plon=-120; Plat=35; CN="LA"
plon=-70.2553 ;plat=43.6615; Plon=-70; Plat=45; CN="Portland,ME"
plon=-122.3321 ;plat=47.6062; Plon=-120; Plat=50; CN="Seattle"
plon=-122.9007 ;plat=47.0379; Plon=-125; Plat=45; CN="Olympia"
plon=-84.5555 ;plat=42.7325; Plon=-80; Plat=45; CN="Lansing"



# plon=-77.0369
# plat=38.9072
idLON=(plon-2*5)<=LON & LON<(plon+2*5)
idLAT=(plat-2*5)<=LAT & LAT<(plat+2*5)
idlon=(plon-2*1/24)<=lon & lon<(plon+2*1/24)
idlat=(plat-2*1/24)<=lat & lat<(plat+2*1/24)

id0=373:(372+360) # training index
id1=(372+360+1):828 # testing index

x0=X[idLON,idLAT,id0]
x1=X[idLON,idLAT,id1]
y0=y_train[idlon,idlat,]
y1=y_test[idlon,idlat,]

px=sum(idLON)*sum(idLAT)
py=sum(idlon)*sum(idlat)
n0=dim(y0)[3] # training sample size
n1=dim(y1)[3] # testing sample size

X0=matrix(0,n0,px)
X1=matrix(0,n1,px)
Y0=matrix(0,n0,py)
Y1=matrix(0,n1,py)

for(i in 1:n0){
  X0[i,]=c(x0[,,i])
  Y0[i,]=c(y0[,,i])
}

for(i in 1:n1){
  X1[i,]=c(x1[,,i])
  Y1[i,]=c(y1[,,i])
}

X0=cbind(rep(1,n0),X0)
X1=cbind(rep(1,n1),X1)
px=px+1
#----------------------
# normalize <- function(newdataf, dataf){
#   normalizeddataf <- newdataf 
#   for (n in 1:ncol(newdataf)){
#     normalizeddataf[,n] <-  
#       (newdataf[,n] - mean(dataf[,n])) /  (sd(dataf[,n]) )
#   } 
#   return(normalizeddataf)
# }
# 
# Y01 <- normalize(Y0, Y0) 
# Y11  <- normalize(Y1, Y0) 
# 
# Y0=Y01
# Y1=Y11
#----------------------

X0l=X1l=list()
for(i in 1:py){
  X0l[[i]]=X0
  X1l[[i]]=X1
}

data=list(Xt=as.matrix(bdiag(X0l)), y=c(Y0)) 

ctr=list(L = 1, use.gram = TRUE, maxiter = 5000, tol = 1e-2, init = NULL, sigma = 0.9)

#Lpx=crossprod(fd(px)); Lpy=crossprod(fd(py))
Lpx=matrix(0,px,px); Lpx[-1,-1]=gridlap(4,4); Lpy=gridlap(4,4)
#Lpx=matrix(0,px,px); Lpx[-1,-1]=cor(X0[,-1]); Lpy=gridlap(4,4)

#Lpx=cor(X0);Lpx[is.na(Lpx)]=0;Lpy=cor(Y0)
Lx=kronecker(Lpx,diag(rep(1,py))); Ly=kronecker(Lpy,diag(rep(1,px)))
P=matrix(0,px*py,px*py)
for(i in 1:px){ for(j in 1:py){ P[(i-1)*py+j,(j-1)*px+i]=1 } } # Permutation matrix

lambda1=0.006602346
lambda2=0
lambdaly=0.5
Lambda=lambdaly*Ly+(1-lambdaly)*t(P)%*%Lx%*%P
Lambda=diag(rep(1,ncol(data$Xt)))
stepsize="fixed" #"backtracking" #
stepsize="backtracking"; ctr$L=NULL
#ctr$init=rnorm(ncol(data$Xt),0,10)
ctr$sigma=0.9
#ctr$L=ifelse(stepsize=="fixed",Ll2(data$Xt, Lambda, lambda2=lambda2),NULL)

out = senet_fista(data$Xt, data$y, Lambda, lossfun=l2loss, gradfun=l2grad, lambda1=lambda1, lambda2=lambda2, stepsize=stepsize, control = ctr)
Bhat=matrix(out$beta,byrow = FALSE,nrow = px, ncol = py)

# t0=proc.time()
# outcv=cv.msenet_fista(nfold=5,py=py,X=data$Xt,y=data$y,Lambda=Lambda,lossfun=l2loss,gradfun=l2grad,lambda1s=c(1,1000),lambda2s=c(1,1000),stepsize="fixed",control=ctr)
# proc.time()-t0
# Bhat=matrix(outcv$beta,byrow = FALSE,nrow = px, ncol = py)

data1=list(Xt=as.matrix(bdiag(X1l)), y=c(Y1)) 

Y_test_hat=as.matrix(bdiag(X1l))%*%(out$beta)

E_test=sqrt(mean((c(Y1)-Y_test_hat)^2))
E_test
##### other solver ####
cvfit2=cv.glmnet(data$Xt, data$y,nfolds = 5,alpha=0)
Y1hat2=predict(cvfit2,as.matrix(bdiag(X1l)))
sqrt(mean((Y1hat2-c(Y1))^2))

# #library(quadrupen)
# Leig <- eigen(Lambda); V <- Leig$vectors; lam <- Leig$values
# if(any(lam<0)){
#   Lambda=V%*%diag(lam*(lam>=0))%*%solve(V)           #,struct=Lambda+(1e-5)*diag(rep(1,ncol(data$Xt)))
# }
# 
# fitquad=quadrupen::elastic.net(x=data$Xt,y=data$y,lambda1=0.0006,lambda2=0,normalize = T,control=list(method="fista",checkargs=F))
# Bhatquad=fitquad@coefficients
# sqrt(mean((as.matrix(bdiag(X1l))%*%t(Bhatquad)-c(Y1))^2))
#######################

# LASSO
Elasso=0
for(i in 1:py){
  cvfit=cv.glmnet(X0, Y0[,i],nfolds = 5,alpha=1)
  Y1hat=predict(cvfit,X1)
  Elasso=Elasso+sum((Y1hat-Y1[,i])^2)
}
Elasso=sqrt(Elasso/(n1*py))
Elasso

########### BCSD ###########
#load("../clinet/data/air_mon_mean_mon_mean_not_removed_sub.RData")
#load("../clinet/data/air_mon_mean_mon_mean_removed_sub.RData")

bx0=X[LON==Plon,LAT==Plat,id0]+273.15
bx1=X[LON==Plon,LAT==Plat,id1]+273.15
by0=y_train[idlon,idlat,]
by1=y_test[idlon,idlat,]

Yim=apply(by0,3,mean)
Yim_test=apply(by1,3,mean)

Xim=bx0
Xim_test=bx1

Qi=quantile(Yim,ecdf(Xim)(Xim))
Qi_test=quantile(Yim,ecdf(Xim)(Xim_test))

# plot(sort(Xim),rank(sort(Xim))/length(Xim),type="o",xlim=c(min(c(Yim,Xim)),max(c(Yim,Xim))),xlab="T",ylab="ECDF",
#      main=paste("Quantile Mapping (",CN,")"),col="blue")
# lines(sort(Yim),rank(sort(Yim))/length(Yim),type="o",col="red")
# legend("bottomright",c("GCM","Obs"),lty=1,pch=1,col=c("blue","red"))

rmse2=rep(0,96)

for(k in 1:96){
  #T1=apply(Pr1_test[,,(1:Mon[mon,3])+(k-1)*Mon[mon,3]],c(1,2),sum, na.rm = TRUE) # Test truth (not available in real case) (for a given year)
  T1=by1[,,k]
  #M=apply(Pr1,c(1,2),sum, na.rm = TRUE)/25 # Local Historical Monthly (accumulation) mean
  M=apply(by0,c(1,2),mean)
 D2=Qi_test[k]*M/mean(Qi) #QM Downscale 1979 Jan, actually apply to training data (which need to try in testing data)
 rmse2[k]=sqrt(mean((D2-T1)^2)) 

}

E_BCSD=sqrt(mean(rmse2^2))#sqrt(sum(rmse2^2))
E_BCSD
#=======================
TES=0
for(k in 1:96){
  T1=by1[,,k]
  M=apply(by1,c(1,2),mean)
  D2=Qi_test[k]*M/mean(Qi)
  rmse2[k]=sqrt(mean((D2-T1)^2)) 
  TES=TES+sum((D2-T1)^2)
  
}

sqrt(TES/(n1*py))

# plot(1:96,rmse2,col="black",type="o",ylim = c(0,max(c(rmse2,rmse3))),xlab="",ylab="RMSE",
#      main=paste("Downscaling on Testing Set in",Mon[mon,1]))
# lines(1:96,rmse3,col="green",type="o")



#df=NULL

df=rbind(df,c(CN,E_test,Elasso,E_BCSD))

colnames(df)=c("City","DSSMR","LASSO","BCSD")
save(df,file="results/real_data_l1000l1_alluse_removed.RData")

# ########## Real Data plot ############
# 
# load("results/real_data_l10.RData")
# dfr=data.frame(x=c(df[,"City"],df[,"City"]),y=round(as.numeric(c(df[,"DSSMR"],df[,"BCSD"])),2),type=as.factor(c(rep(" DSSMR",6),rep("BCSD",6))),od=c(rep(2,6),rep(1,6)))
# ggplot(dfr, aes(x=x, y=y, fill=type)) +
#   geom_bar(stat="identity", position=position_dodge())+
#   xlab("City")+ylab("RMSE")+ggtitle("Downscaling Results")+theme(plot.title = element_text(hjust = 0.5))
# 
# ########## Real Data map plot ############
# library(ggmap)
# library(ggplot2)
# 
# plon=-71.0589 ;plat=42.3601; Plon=-70; Plat=40 ;CN="Boston"
# plon=-71.4676 ;plat=42.7654; Plon=-70; Plat=45 ;CN="Nashua"
# plon=-71.4128 ;plat=41.8240; Plon=-70; Plat=40 ;CN="Providence"
# plon=-74.0060 ;plat=40.7128; Plon=-75; Plat=40 ;CN="New York"
# plon=-77.0369 ;plat=38.9072; Plon=-75; Plat=40 ;CN="Washington, D.C."
# plon=-72.6851 ;plat=41.7637; Plon=-75; Plat=40 ;CN="Hartford"
# 
# d=data.frame(lat=c(42.3601,42.7654,41.8240,40.7128,38.9072,41.7637),
#              lon=c(-71.0589,-71.4676,-71.4128,-74.0060,-77.0369,-72.6851),
#              city=c("Boston","Nashua","Providence","New York","Washington, D.C.","Hartford"))
# 
# d <- data.frame(lat=c(50.659631, 50.607213, 50.608129),
#                 lon=c(3.09319, 3.011473, 3.031529))
# 
# Boston <- get_map("New York", zoom=7)
# p <- ggmap(Boston)
# p + geom_point(data=d, aes(x=lon, y=lat), color="red", size=10, alpha=0.5)
# 







