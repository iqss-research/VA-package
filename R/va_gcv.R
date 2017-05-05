## a cross-validation method for finding best number of symptoms

va.gcv<-function(formula, data, nsymp.vec, n.subset=300, prob.wt=1,
  boot.se=FALSE, nboot=1, printit=FALSE, print.reg.size=TRUE)
  {
    mse<-rep(0, length(nsymp.vec))
       for (i in 1:length(nsymp.vec)) {
         res<-va(formula, data, nsymp=nsymp.vec[i],
                      n.subset=n.subset, boot.se=boot.se,
                      nboot=nboot, printit=printit,
                      prob.wt=prob.wt,
                      print.reg.size=print.reg.size)
         mse[i]<-mean((res$true.CSMF-res$est.CSMF)^2)
       }                      

    return(list(best.nsymp=nsymp.vec[mse==min(mse)], mse=mse))



  }
