############ Default parameters ############
ctr=list(lambda1=1000, lambda2=1000,muy=1, stepsize="fixed",lx="1d",ly="1d", 
         L = 1, use.gram = TRUE, maxiter = 2000, tol = 1e-5, init = NULL, sigma = 0.9)
c1=list(n=200,n_test=100,px=40,py=20,noise.signal=0.5,noise.response=sqrt(5),dataseed=1,typeb="onebump",typeB="circ",typeX="5signal")

############ Change parameters ############
#l1s=c(0,1,10,100,1000,10000)
#l2s=c(0,1,10,100,1000,10000)
l1s=1000
l2s=1000
lys=c(0,0.5,1) #(0:10)/10 #
ns= (1:5)*100 #c(200,400,600) 

df=NULL
for(i1 in 1:length(l1s)){
  for(i2 in 1:length(l2s)){
      for(i3 in 1:length(lys)){
          for(i4 in 1:length(ns)){
               df=rbind(df,c(i1,i2,i3,i4,l1s[i1],l2s[i2],lys[i3],ns[i4]))
          }
      }
  }
}
colnames(df)=c("id1","id2","id3","id4","lambda1","lambda2","lambdaly","n")

############ Simulation function ############ 
simu=function(jobid,ctr,c1,df){
    # 1. Update parameters
    ctr$muy=df[jobid,"lambdaly"]
    ctr$lambda1=df[jobid,"lambda1"]
    ctr$lambda2=df[jobid,"lambda2"]
    c1$n=df[jobid,"n"]
    # 2. Generate Data
    data=getDataSim(n=c1$n,n_test=c1$n_test,px=c1$px,py=c1$py,
                    noise.signal=c1$noise.signal,noise.response=c1$noise.response,dataseed=c1$dataseed,
                   typeb=c1$typeb,typeB=c1$typeB,typeX=c1$typeX)
    # 3. Run algo
    res=msnet_fista(data,ctr=ctr,jobid=jobid) 
}