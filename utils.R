simblock1d <- function(n, p=100, noise.signal = 0.25, noise.response = 30, beta,dataseed=1,...){
  set.seed(dataseed)
  Xt <- matrix(ncol = p, nrow = n)
  for(i in 1:n){
    bi <- runif(5, 0, 5)
    mi <- runif(5, 0, 2*pi)
    Xt[i,]  <- sapply(1:p, function(z) sum(bi * sin(z*pi*(5-bi)/50-mi))) + rnorm(p, sd=noise.signal)
  }
  y <- Xt %*% beta + rnorm(n, sd=noise.response)
  return(list(Xt = Xt, y=y, beta = beta))
}
                      
                      
fd <- function(p){
  D <- matrix(nrow = p-1, ncol = p, data = 0)
  for(i in 1:(p-1)){
    D[i,i] <- 1
    D[i,i+1] <- -1
  }
  D
} 

simblock1d2 <- function(n, p=100, noise.signal = 0.25, noise.response = 30, beta,dataseed=1,...){
  set.seed(dataseed)
  Xt <- matrix(ncol = p, nrow = n)
  for(i in 1:n){
    Mp=rep(0,p)  
    Ip=diag(rep(1,p))
    Xt[i,] = mvrnorm(1,Mp,Ip)
  }
  y <- Xt %*% beta + rnorm(n, sd=noise.response)
  return(list(Xt = Xt, y=y, beta = beta))
}                      
                      
                      
twobump=function(t,px){
  (floor(0.2*px) < t && t < floor(0.4*px))*((-1)*(floor(0.3*px) - t)^2+px)/(2*px) +
  (floor(0.6*px) < t && t < floor(0.8*px))*((floor(0.7*px) - t)^2-px)/(2*px)
}

onebump=function(t,px){
  (floor(0.2*px) < t && t < floor(0.4*px))*((-1)*(floor(0.3*px) - t)^2+px)/(2*px)
  (floor(0.2*px) < t && t < floor(0.4*px))*(1-abs((floor(0.3*px) - t)/floor(0.1*px)))  
}
# a=sapply(1:px, function(t) (20 < t && t < 40)*((-1)*(30 - t)^2+100)/200 +
#          (60 < t && t < 80)*((70 - t)^2-100)/200)

fourblock=function(t,px){
  (floor(0.2*px) < t && t <= floor(0.3*px))*0.5+
  (floor(0.3*px) < t && t <= floor(0.4*px))*1+
  (floor(0.4*px) < t && t <= floor(0.5*px))*0.5+
  (floor(0.5*px) < t && t <= floor(0.6*px))*0.25
}

# Generate 2d laplacian matrix
genlap=function(p1,p2){
  n=p1*p2
  L=W=D=matrix(0,n,n)
  
  for(i1 in 1:p1){
    for(j1 in 1:p2){
      for(i2 in 1:p1){
        for(j2 in 1:p2){

          x=abs(i1-i2)
          y=abs(j1-j2)          
          
          #logi=((x==1)&(y==0))|((x==0)&(y==1))
          logi=((x+y)==1)
            
          if(logi){
            i=(i1-1)*p2+j1 #label by row
            j=(i2-1)*p2+j2
            W[i,j]=1
          }
        }
      }
    }
  }
  M=apply(W,1,sum) 
  D=diag(M)
  L=D-W
  return(list(L=L,D=D,W=W))
}
                      
                      
# Generate Synthetic data
getDataSim=function(n,n_test,px,py,noise.signal, noise.response,dataseed=1,typeb="twobump",typeB="circ",typeX="5signal"){
    # generate simulated data (X,B,Y) for msnet
    # typeb: how to generate b - bump, block
    # typeB: how to generate B - delay(circle), scale
    # typeX: how to generate X - in paper each feature is combo of 5 signals, each feature is independent
    #============= 1. B (coefficients) =============
    # 1.1 generate b
    if(typeb=="twobump"){
        betastar=sapply(1:px,twobump,px)
    }else if(typeb=="onebump"){
        betastar=sapply(1:px,onebump,px)
    }else if(typeb=="fourblock"){
        betastar=sapply(1:px,fourblock,px)
    }
    # 1.2 generate B (according to b)
    B=matrix(0,px,py) 
    if(typeB=="circ"){
        for(i in 1:py){
            B[,i]=betastar[(1:px-i)%%px+1]
        }
    }else if(typeB=="quadscale"){
        for(i in 1:py){
            B[,i]=betastar*2*(0.5+abs(i-py/2)/py) 
        }
    }else if(typeB=="quadscale2"){
        for(i in 1:py){
            if(i<py/2){B[,i]=betastar*(0.5+abs(i-py/2)/py)}else{B[,i]=-betastar*(0.5+abs(i-py/2)/py)}
            #B[round(py/2):py,]=0
        }
    }    
    # now we have B: px*py
    #============= 2. X (and Y) =============
    Y=matrix(0,n,py); X=list(py); Y_test=matrix(0,n_test,py); X_test=list(py);     
    if(typeX=="5signal"){
        for(i in 1:py){
            # trainging
#            data_wide=simblock1d(n=n, p=px, noise.signal = noise.signal, noise.response = noise.response,beta=B[,i],dataseed=dataseed)
#            Y[,i]=data_wide$y; X[[i]]=data_wide$Xt
            # testing
#            data_wide_test=simblock1d(n=n_test, p=px, noise.signal = noise.signal, noise.response = noise.response, beta=B[,i],dataseed=dataseed+20170223)
#            Y_test[,i]=data_wide_test$y; X_test[[i]]=data_wide_test$Xt
            
            data_wide_all=simblock1d(n=n+n_test, p=px, noise.signal = noise.signal, noise.response = noise.response,beta=B[,i],dataseed=dataseed)
            Y[,i]=data_wide_all$y[1:n]; X[[i]]=data_wide_all$Xt[1:n,]
            Y_test[,i]=data_wide_all$y[(n+1):(n+n_test)]; X_test[[i]]=data_wide_all$Xt[(n+1):(n+n_test),]
            
        }
    }else if(typeX=="ind"){
        for(i in 1:py){
            # trainging
            data_wide=simblock1d2(n=n, p=px, noise.signal = noise.signal, noise.response = noise.response,beta=B[,i],dataseed=dataseed)
            Y[,i]=data_wide$y; X[[i]]=data_wide$Xt
            # testing
            data_wide_test=simblock1d2(n=n_test, p=px, noise.signal = noise.signal, noise.response = noise.response, beta=B[,i],dataseed=dataseed+20170223)
            Y_test[,i]=data_wide_test$y; X_test[[i]]=data_wide_test$Xt
        }
    }
    # now have: B, X, Y, X_test, Y_test 
    #============= 3. save B, X, Y =============    
    data=list(B=B,X=X,Y=Y,X_test=X_test,Y_test=Y_test)    
    return(data)
}
    
# Main algo    
msnet_fista=function(data,ctr=ctr,jobid=0,df){
    # take input data, and do multi-task structured enet
    #============= 1. Wide data to long data =============
    # data=list(B=B,X=X,Y=Y,X_test=X_test,Y_test=Y_test) 
    # c(A): stack column of A together
    Data=list(X=as.matrix(bdiag(data$X)), y=c(data$Y), beta = c(data$B),
              X_test=as.matrix(bdiag(data$X_test)),y_test=c(data$Y_test),
              p1x=data$p1x,p2x=data$p2x,p1y=data$p1y,p2y=data$p2y) # same X # npy*pxpy
    px=dim(data$B)[1]; py=dim(data$B)[2]; n=dim(data$X[[1]])[1]; n_test=dim(data$X_test[[1]])[1]; 
    #============= 2. Quadratic Matrices Lx, Ly =============
    lx=ctr$lx; ly=ctr$ly; muy=ctr$muy;
    #------------- 2.1 Lpx, Lpy wide version-------------
    # Lpx
    if(lx=="1d"){
        Lpx=crossprod(fd(px))
    }else if(lx=="2d"){
        Lpx=genlap(p1x,p2x)$L #p1: num of lat, p2: num of lon
    }else if(lx=="cor"){
        Lpx=cor(data$X[[1]])
    }
    # Lpy
    if(lx=="1d"){
        Lpy=crossprod(fd(py))
    }else if(lx=="2d"){
        Lpy=genlap(p1y,p2y)$L #p1: num of lat, p2: num of lon
    }else if(lx=="cor"){
        Lpy=cor(data$Y)
    }    
    #------------- 2.2 Lx(transform), Ly long version-------------
    Lx=kronecker(Lpx,diag(rep(1,py)))
    Ly=kronecker(Lpy,diag(rep(1,px)))
    P=matrix(0,px*py,px*py)
    for(i in 1:px){ for(j in 1:py){ P[(i-1)*py+j,(j-1)*px+i]=1 } } # Permutation matrix
    Lx=t(P)%*%Lx%*%P
    #------------- 2.3 Combine Lx and Ly-------------
    Lambda=muy*Ly+(1-muy)*Lx
    #============= 3. Fit the model, using senet =============
    lambda1=ctr$lambda1
    lambda2=ctr$lambda2
    stepsize=ctr$stepsize
    ctr$L=ifelse(stepsize=="fixed",Ll2(Data$X, Lambda, lambda2=lambda2),NULL)
    #ctr$L=2*(ctr$L)
    out = senet_fista(Data$X, Data$y, Lambda, lossfun=l2loss, gradfun=l2grad, 
                      lambda1=lambda1, lambda2=lambda2, stepsize=stepsize, control = ctr) 
    Bhat=matrix(out$beta,byrow = FALSE,nrow = px, ncol = py)    
    #============= 4. Prediction using fitted model =============
    Y_test_hat=as.matrix(bdiag(Data$X_test))%*%(out$beta)
    E_test=sqrt(mean((c(Data$y_test)-Y_test_hat)^2))    
    #============= 5. Save the results =============
    slist=list(fit=out,input=list(B=data$B,px=px,py=py,n=n,Lambda=Lambda,lambda1=lambda1,lambda2=lambda2,muy=muy),output=list(Bhat=Bhat,rmse_test=E_test,Y_test=Data$y_test,Y_test_hat=Y_test_hat,obj=out$obj))
    fname=paste0("./results/n",n,"px",px,"py",py,"job",jobid)
    save(slist,file = paste0(fname,".RData"))   
    return()    
}    
    
    
                      