data(VAdata)

all.true<-table(VAdata$cod)

 DSS<-VAdata[VAdata$region==1,-1]
 FBA<-VAdata[VAdata$region==2,-1]

 index.all<-names(all.true)
  index.all<- index.all[order(as.numeric(index.all))]
  comp<-matrix(0, length(index.all),2)
  comp<-as.data.frame(comp)
  dimnames(comp)<-list(index.all, c("true", "emp"))


 res1<-va(formula=cbind(S1+...+S49)~cod,
          data=list(FBA,DSS),  nsymp=16,
                n.subset=5,printit=TRUE,boot.se=FALSE)

  cod.emp<-res1$est.CSMF

  ##compare results
  ##true value
  cod.true<-(res1$true.CSMF)

   plot(cod.emp, cod.true,cex=0.8, xlab="True", ylab="Estimate", main="cause-specific mortality estimation", xlim=c(0,0.3), 
ylim=c(0, 0.3), pch=19, col="red")
  abline(0,1)











 
