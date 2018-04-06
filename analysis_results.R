########## plot what are we doing ############
library(fields)
load("results/n2-1000pxpy1-40wlasso/n2000px40py40seed50job1000Tb.RData")
par(mfrow=c(1,1))
plot(slist$input$B[1,],col="blue",xlab="feature index",ylab = "beta",main="coefficients in one task")
image.plot(slist$input$B, main="True coefficients",xlab="y",ylab="x")
image.plot(slist$output$Bhat,main="Estimation with DSSMR")
image.plot(slist$output$Bhatlasso,main="Estimation with Lasso")

colorbar.plot(1,0.5, 1:10,horizontal = F,col=rgb(1,(1:10)/10,(1:10)/10))

image(as(slist$input$B, "dgCMatrix"))
image(as(slist$output$Bhat, "dgCMatrix"))
image(as(slist$output$Bhatlasso, "dgCMatrix"))

########## simulation table & plot results ############
library(ggplot2)
load("./results/n2-1000pxpy1-40wlasso/df.RData")
load("./results/n2-1000pxpy1-40wlassolambdalasso/df.RData")
#df[,"ferr"]<df[,"ferrlasso"]

df0=as.data.frame(df[,c("n","p","ferr","ferrlasso","dataseeds")])
df0m=aggregate(cbind(df0$ferr,df0$ferrlasso),by=list(n=df0$n,p=df0$p),FUN=mean)
#aggregate(cbind(df0$ferr,df0$ferrlasso),by=list(n=df0$n,p=df0$p),FUN=sd)
colnames(df0m)=c("n","p","ferr","ferrlasso")

# plot sample_size vs error
df0=subset(df0m,p==10)
df00=data.frame(x=c(df0[,"n"],df0[,"n"]),y=c(df0[,"ferr"],df0[,"ferrlasso"]),type=c(rep("DSSMR",5),rep("Lasso",5)))
ggplot(df00, aes(x=x, y=y, fill=type)) + 
  geom_bar(stat="identity", position=position_dodge())+ 
  xlab("Sample Size")+ylab("RMSE")+ggtitle("Estimation Results")+theme(plot.title = element_text(hjust = 0.5))

# plot dimension vs error
df0=subset(df0m,n==2000)
df00=data.frame(x=c(df0[,"p"],df0[,"p"]),y=c(df0[,"ferr"],df0[,"ferrlasso"]),type=c(rep("DSSMR",4),rep("Lasso",4)))
ggplot(df00, aes(x=x, y=y, fill=type)) + 
  geom_bar(stat="identity", position=position_dodge())+ 
  xlab("Dimension")+ylab("RMSE")+ggtitle("Estimation Results")+theme(plot.title = element_text(hjust = 0.5))


# +
#   annotate("text", x=3, y=250, label= "Overall MSE (Test)",fontface =2)+
#   annotate("text", x=3, y=230, label= "Double Taper: 82.46\n Sample Cov: 117.06",fontface =1)

# change to wide table
#dft=subset(df,select = c("n","p","ferr","ferrlasso"))
dftw=NULL
for(i in 1:4){
  dftw=cbind(dftw,as.matrix(df0m)[(i*5-4):(i*5),3:4])
}
round(dftw,2)
write.table(round(dftw,2),file="figs/dftw.csv",row.names = F)






