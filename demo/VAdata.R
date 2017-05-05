data(VAdata0)

#all.true<-table(VAdata0$cod)

 DSS<-VAdata0[VAdata0$region==1,-1]
 FBA<-VAdata0[VAdata0$region==2,-1]

# index.all<-names(all.true)
 # index.all<- index.all[order(as.numeric(index.all))]
  #comp<-matrix(0, length(index.all),2)
  #comp<-as.data.frame(comp)
  #dimnames(comp)<-list(index.all, c("true", "emp"))

 res1<-va(formula=cbind(S1+...+S49)~cod,
          data=list(FBA,DSS),  nsymp=16,
                n.subset=10,prob.wt=1,printit=TRUE,boot.se=FALSE)

  cod.emp<-res1$est.CSMF

  ##compare results
  ##true value
  cod.true<-(res1$true.CSMF)

   plot(cod.emp, cod.true,cex=0.8, xlab="True", ylab="Estimate", main="cause-specific mortality estimation", xlim=c(0,0.3), 
ylim=c(0, 0.3), pch=19, col="red")
  abline(0,1)

  res2 <- va(formula=cbind(S1+...+S49)~cod, fix=c("a13=0.216",
          "a18=0.08"),bound=c("0.005<a14<0.01"),
          data=list(FBA,DSS),  nsymp=16,
                n.subset=10,prob.wt=1,printit=TRUE,boot.se=FALSE)






 
