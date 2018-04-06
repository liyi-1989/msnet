load("./results.df.RData")
px=20
py=20
type="b"
ferr=ferrlasso=ferr2=ferr2lasso=TP=FP=TN=FN=TPlasso=FPlasso=TNlasso=FNlasso=FS=FSlasso=rep(0,nrow(df))
for(jobid in 1:nrow(df)){
  n=df[jobid,"n"]
  fname=paste0("./results/n2-1000pxpy1-40wlasso/n",n,"px",df[jobid,"p"],"py",df[jobid,"p"],"seed",df[jobid,"dataseeds"],"job",jobid,"T",type)
  
  if(jobid%%4!=0){
    load(paste0(fname,".RData"))
  }
  
  
  
  ferr[jobid]=base::norm(slist$input$B-slist$output$Bhat,"F")
  ferrlasso[jobid]=base::norm(slist$input$B-slist$output$Bhatlasso,"F")
  
  ferr2[jobid]=base::norm(slist$input$B-slist$output$Bhat,"2")
  ferr2lasso[jobid]=base::norm(slist$input$B-slist$output$Bhatlasso,"2")
  
  TP[jobid]=sum(slist$input$B!=0 & slist$output$Bhat!=0)
  FP[jobid]=sum(slist$input$B==0 & slist$output$Bhat!=0)
  TN[jobid]=sum(slist$input$B==0 & slist$output$Bhat==0)
  FN[jobid]=sum(slist$input$B!=0 & slist$output$Bhat==0)
  FS[jobid]=2*TP[jobid]/(2*TP[jobid]+FN[jobid]+FP[jobid])
  
  TPlasso[jobid]=sum(slist$input$B!=0 & slist$output$Bhatlasso!=0)
  FPlasso[jobid]=sum(slist$input$B==0 & slist$output$Bhatlasso!=0)
  TNlasso[jobid]=sum(slist$input$B==0 & slist$output$Bhatlasso==0)
  FNlasso[jobid]=sum(slist$input$B!=0 & slist$output$Bhatlasso==0)
  FSlasso[jobid]=2*TPlasso[jobid]/(2*TPlasso[jobid]+FNlasso[jobid]+FPlasso[jobid])
}
df=cbind(df,ferr,ferrlasso,ferr2,ferr2lasso,TP,FP,TN,FN,TPlasso,FPlasso,TNlasso,FNlasso,FS,FSlasso)
save(df,file="./results/df.RData")


