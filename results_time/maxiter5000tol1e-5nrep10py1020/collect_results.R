library(plyr)
library(ggplot2)

load("./results_time/maxiter5000tol1e-5nrep10/results_df.RData")

T1=rep(0,nrow(df))
for(i in 1:nrow(df)){
  n=df[i,"n"]; px=df[i,"px"]; py=df[i,"py"]; dataseed=1
  fname=paste0("./results_time/maxiter5000tol1e-5nrep10/n",n,"px",px,"py",py,"seed",dataseed,"job",i)
  load(paste0(fname,".RData"))
  T1[i]=slist$time
}

#df=cbind(df,1:nrow(df))
df=cbind(df,T1)
colnames(df)[9]="time"
df=as.data.frame(df)

#dfv1=aggregate(time ~ py+px+n, as.matrix(df[,5:9]), mean)
df[,"time"]=df[,"time"]/60
dfa=ddply(df,.(n,px,py),summarise,V1=mean(time),V2=sd(time))

#======= py=10

dfplot=subset(dfa,py==10,select = c("n","px","V1","V2"))
dfplot[,"n"]=as.factor(dfplot[,"n"])

ggplot(dfplot, aes(x=px, y=V1, colour=n)) +
  geom_errorbar(aes(ymin=V1-V2, ymax=V1+V2),width=1) +
  geom_line() +
  geom_point()+
  xlab("px") +
  ylab("Time (min)") +
  ggtitle("Time Comparison")+
  theme(plot.title = element_text(hjust = 0.5))

#======= py=20

dfplot=subset(dfa,py==20,select = c("n","px","V1","V2"))
dfplot[,"n"]=as.factor(dfplot[,"n"])

ggplot(dfplot, aes(x=px, y=V1, colour=n)) +
  geom_errorbar(aes(ymin=V1-V2, ymax=V1+V2),width=1) +
  geom_line() +
  geom_point()+
  xlab("px") +
  ylab("Time (min)") +
  ggtitle("Time Comparison")+
  theme(plot.title = element_text(hjust = 0.5))