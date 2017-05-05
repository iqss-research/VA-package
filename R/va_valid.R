## this program implements the symptom select method
## discussed in King and Lu (2008b) WHO report

require(Design, quietly=TRUE)
require(quadprog, quietly=TRUE)

va.validate<-function(formula,  data=list(hospital=NA, community=NA),
               nsymp=10,n.subset=300, nboot=1, boot.se=FALSE, method="quadOpt",
                      fix=NA, bound=NA,
                      prob.wt=1,
                      printit=TRUE,print.reg.size=TRUE,
                      clean.method="ttest",
                      min.Symp=10, confidence=0.95, FDR=TRUE){

if (min.Symp<nsymp) min.Symp<-nsymp

if (method=="quadOpt") {
    tt<-requireNamespace("quadprog")
    if (!tt)
      stop("\n Require Quadractic Programming package 'quadprog'\n",
           call.=FALSE)
  }

  prior <- list()

  if (!is.na(fix)) {
    prior$fix <- fFix(fix)
  }
  if (!is.na(bound)) {

    prior$bdd <- fBdd(bound)
  }
  print(prior)

  rhsVars<-all.vars(formula[[3]])
  lhsVars<-all.vars(formula[[2]])
  lhsVars<-dot(lhsVars,data[[1]])

 if(length(rhsVars)==2){
    hospVars<-c(rhsVars[[1]],lhsVars)
    commVars<-c(rhsVars[[2]],lhsVars)
  }
  else
    hospVars<-commVars<-c(rhsVars,lhsVars)

  ## first column cause of death
  ## second to last symptoms

  hospD <- data[[1]][,hospVars]


 comm.cod <- TRUE
  commD <- data[[2]][,commVars]
  if ((rhsVars[[1]] %in% colnames(commD))==FALSE){
    warning("cause of deach not in population sample")
    comm.cod <- FALSE
    tmp <- rep(1,nrow(commD))
    commD<-cbind(tmp,commD)
    colnames(commD)[[1]]<-rhsVars[[1]]
  }

  q<-ncol(hospD)
  p<-nrow(hospD)

  ##check data: exit if
  ## any missing in causes of death (column 1)
  ## or, any missing in symptoms (column 2:q)
  ## or, no variation in certain symptoms (column 2:q)
  ## or, any values not equal to 0 or 1 (column 2:q)

  if (any(is.na(hospD[,1])))
    stop("\n Missing values found in the cause of death column\n",
         call.=FALSE)

  if (any(is.na(hospD[,2:q])) || any(is.na(commD[,2:q])))
    stop("\n Missing valus found in symptom profiles\n", call.=FALSE)

  colsums<-colSums(hospD[,2:q], na.rm=TRUE)

  if (any(colsums==0 | colsums ==p))
    stop("\nSome symptoms have no variation!!\n",call.=FALSE)

  if ((!all(as.matrix(hospD[,2:q]) %in% c(0,1))) ||
      (!all(as.matrix(commD[,2:q]) %in% c(0,1))))
    stop("\n All te symptoms should be coded 0,1 \n",call.=FALSE)


  ##if user decides to supply prob.wt, it must be the same length as
  ##symptoms

  if (length(prob.wt)>1)
    if (length(prob.wt)!=dim(hospD)[2]-1)
      stop("probability weights
  length must equal to number of symptoms in hosptial/community
          data\n")


  cod.names<-union(names(table(hospD[,1])), names(table(commD[,1])))

  index.all<-cod.names
  index.all<- index.all[order(as.numeric(index.all))]
  nd <- length(index.all)
  ns <- ncol(hospD)-1

  if (length(prob.wt)==1) {
    if (prob.wt==1) {
      temp<-rbind(hospD, commD)
      end<-dim(temp)[2]
      sta<-2
      prob.wt<-rep(0, (end-sta+1))
      for ( i in sta:end){
        tt<-table(temp[,i])
        n<-sum(tt)
        prob.wt[i-sta+1]<-tt[1]*tt[2]/n^2
      }
      prob.wt<-prob.wt/sum(prob.wt)
    } else {
      if (prob.wt==0) prob.wt<-NULL
    }
  }


    res<- cod.est.base(hospital=hospD,community=commD,
                       n.subset=n.subset, nsymp=nsymp,
                       method=method,prior=prior,
                       prob.wt=prob.wt,printit=printit,
                       print.reg.size=print.reg.size,
                       checkrow.names=FALSE, logfile="junk.txt")

cod.emp<-cod.true<-rep(0,nd)
    cod.emp[index.all %in% names(res$est.CSMF)]<-res$est.CSMF
    cod.true[index.all %in% names(res$true.CSMF)]<-res$true.CSMF

    ##P(S|D) in hospital
    sprob.h<-matrix(0, nd, ns)
    for (i in 1:nd)
      sprob.h[i,]<-colMeans(hospD[hospD[,1]==index.all[i],-1])


###calcualte variance in hat{P}(S) based on bootstrap
if (nboot==1) boot.se<-FALSE
if (boot.se) {
      comp<-matrix(0, nboot, nd)

      s.fit <- e.boot<-matrix(0, nboot, ncol(hospD)-1)
      for (t in 1:nboot){
        p1<-dim(commD)[1]
        p2<-dim(hospD)[1]

        boot.good<-FALSE
        while (!boot.good) {

          index.boot1<-sample(1:p1, p1, replace=TRUE)
          index.boot2<-sample(1:p2, p2, replace=TRUE)

          commD0<-commD[index.boot1,]
          hospD0<-hospD[index.boot2,]

          boot.good<-!(any(colMeans(commD0[,-1])==1) ||
                     any(colMeans(commD0[,-1])==0) ||
                     any(colMeans(hospD0[,-1])==1) ||
                     any(colMeans(hospD0[,-1])==0) )

        if (boot.good) {
            res0<- try(cod.est.base(hospital=hospD0,community=commD0,
                       n.subset=n.subset, nsymp=nsymp,
                            method=method,prior=prior,
                       prob.wt=prob.wt,printit=FALSE,
                       print.reg.size=FALSE,
                       checkrow.names=FALSE), TRUE)

        if (length(res0)==1) boot.good<-FALSE }
        }

#        print(res0$est.CSMF)
        comp[t,index.all %in% names(res0$est.CSMF)]<-res0$est.CSMF

        sd.mat<-matrix(0,ns,nd)
        for (k in 1:nd)
          {
            sd.mat[,k] <- colMeans(hospD0[hospD0[,1]==k,2:(ns+1)])
          }
        sd.mat[is.na(sd.mat)]<-0
        s.fit[t,]<-sd.mat%*%comp[t,]
        e.boot[t,]<-s.fit[t,]-colMeans(commD0[,-1])
    if (boot.se)
      cat(paste(t, "th bootstrap sampling done.\n\n"))
      }
    }


    e.var<-hatS.var<-res.out<-res.pred<-delete.list<-res.out.boot<-res.pred.boot<-list()
    res.out[[1]]<-cod.emp
    res.pred[[1]]<-t(sprob.h)%*%cod.emp
    delete.list[[1]]<-NULL
if (boot.se) {
res.out.boot[[1]]<-comp
res.pred.boot[[1]]<-s.fit
     hatS.var[[1]]<-apply(s.fit,2, var)
e.var[[1]]<-apply(e.boot,2,var)
}
    PS<-colMeans(commD[,-1])

    e0<-as.vector(res.pred[[1]]-PS)

cat("clean.method=", clean.method, "\n")
##t test
    if (clean.method=="ttest") {
      ind<-ttest(x=e0, prob=1-(1-confidence)/2, df=ns-1)

    }

if (clean.method=="ztest") {
  ind<-ztest(x=e0, sig2=apply(e.boot,2,var), prob=1-(1-confidence)/2)
  cat("\n ztest\n")
}
## q test
if (clean.method=="qtest") {
      ind<-qtest(e0, prob=1-(1-confidence)/2, size=ns)
          cat("\n qtest\n")
    }

index.S<-NULL
    clean<-TRUE

    if (any(ind)) {
      clean<-FALSE

      j<-2

      out.list<-(1:ns)[abs(e0)==max(abs(e0))]
      index.S<-(1:ns)[abs(e0)!=max(abs(e0))]
      cat("ind \n")
      print(index.S)

    }

    while (!clean) {


        res<-cod.est.base(hospital=hospD[,c(1,(index.S+1))],
                          community=commD[,c(1,(index.S+1))],
                       n.subset=n.subset, nsymp=nsymp,
                       method=method,prior=prior,
                       prob.wt=prob.wt[index.S],printit=printit,
                       print.reg.size=print.reg.size,
                       checkrow.names=FALSE, logfile="junk.txt")
        cod.emp<-rep(0,nd)
        cod.emp[index.all %in% names(res$est.CSMF)]<-res$est.CSMF

        if (boot.se) {
          comp<-matrix(0, nboot, nd)

          s.fit <- e.boot<-matrix(0, nboot, length(index.S))

          for (t in 1:nboot){


            p1<-dim(commD)[1]
            p2<-dim(hospD)[1]


            boot.good<-FALSE
            while (!boot.good) {

              index.boot1<-sample(1:p1, p1, replace=TRUE)
              index.boot2<-sample(1:p2, p2, replace=TRUE)

              commD0<-commD[index.boot1,c(1,(index.S+1))]
              hospD0<-hospD[index.boot2,c(1,(index.S+1))]

              boot.good<-!(any(colMeans(commD0[,-1])==1) ||
                     any(colMeans(commD0[,-1])==0) ||
                     any(colMeans(hospD0[,-1])==1) ||
                     any(colMeans(hospD0[,-1])==0) )

              if (boot.good) {
                res0<- try(cod.est.base(hospital=hospD0,community=commD0,
                       n.subset=n.subset, nsymp=nsymp,
                            method=method,prior=prior,
                       prob.wt=prob.wt[index.S],printit=FALSE,
                       print.reg.size=FALSE,
                       checkrow.names=FALSE, logfile="junk.txt"),TRUE)
              if (length(res0)==1) boot.good<-FALSE }
            }

#            print(res0$est.CSMF)
            comp[t,index.all %in% names(res0$est.CSMF)]<-res0$est.CSMF
            ns0<-ncol(hospD0)-1
            sd.mat<-matrix(0,ns0,nd)
            for (k in 1:nd)
              {
                sd.mat[,k] <- colMeans(hospD0[hospD0[,1]==k,2:(ns0+1)])
              }
            sd.mat[is.na(sd.mat)]<-0
            s.fit[t,]<-sd.mat%*%comp[t,]
            e.boot[t,]<-s.fit[t,]-colMeans(commD0[,-1])
            if (boot.se)
              cat(paste(t, "th bootstrap sampling done.\n\n"))
          }
        }


        res.out[[j]]<-cod.emp
        res.pred[[j]]<-t(sprob.h)%*%cod.emp
        delete.list[[j]]<-out.list

if (boot.se) {
 hatS.var[[j]]<-apply(s.fit,2, var)
        e.var[[j]]<-apply(e.boot,2,var)
        res.out.boot[[j]]<-comp
        res.pred.boot[[j]]<-s.fit
      }
        ee<-as.vector(res.pred[[j]]-PS)[index.S]
       ##t test
        if (clean.method=="ttest") {

          ind<-ttest(x=ee, prob=1-(1-confidence)/2, df=length(ee)-1) }
## ztest

        if (clean.method=="ztest") {
          ind<-ztest(x=ee, sig2=apply(e.boot,2,var), prob=1-(1-confidence)/2)
          cat("\n ztest\n")
}

        ## q test
        if (clean.method=="qtest") {
          ind<-qtest(ee, prob=1-(1-confidence)/2, size=length(ee))
        }


        clean<-TRUE

        if (any(ind)) {
        out.list<-c(out.list, index.S[abs(ee)==max(abs(ee))])
        index.S<-index.S[abs(ee)!=max(abs(ee))]
        clean<-FALSE
        cat("ind \n")
        print(index.S)
      }

        j<-j+1
      if (length(index.S)<=min.Symp) clean<-TRUE
    }

if (FDR) {
  Ps<-colMeans(commD[,-1])
  ns<-length(Ps)
  nr<-length(res.out)
  tscore<-rep(0,nr)
  e0<-as.numeric(res.pred[[1]]-Ps)
tscore[1]<-max(abs(e0-mean(e0))/(var(e0)^0.5))
if (nr>=2) {
  for (k in 2:nr) {
  e0<-(res.pred[[k]]-Ps)[-delete.list[[k]]]
  tscore[k]<-max(abs(e0-mean(e0))/(var(e0)^0.5))
  }
}

tvalue<-qt(1-(1-confidence)/(1:nr), (ns-1):(ns-nr))

}


if (boot.se) {
  out<-list(cod.list=res.out, Ps.list=res.pred, Ps.var=hatS.var, e.var=e.var, delete.list=delete.list, cod.list.boot=res.out.boot, Ps.list.boot=res.pred.boot)
}
else out<-list(cod.list=res.out, Ps.list=res.pred, delete.list=delete.list)


if (FDR & length(delete.list)>=1) {
  k<-1
  while (tscore[k]>tvalue[k]) {
    k<-k+1
    out$FDR.delete.list<-delete.list[[k]]
  }




}
 if (!is.null(index.S) && length(index.S)<=min.Symp) out$list.short<-"yes"

  return(out)
}




ttest<-function(x, prob, df) {
    x<-abs(x-mean(x))/var(x)^0.5
    return(x>qt(prob, df))
  }

ztest<-function(x, sig2, prob) {
x<-abs(x)/(sig2^0.5)
return(x>qnorm(prob))
}


qtest<-function(x, size, prob) {
  qtable<-matrix(0, 6,3)
  qtable[1,]<-c(0.642, 0.710, 0.821)  #size=5
  qtable[2,]<-c(0.412, 0.466, 0.568)  #size=10
  qtable[3,]<-c(0.338, 0.384, 0.475)  #size=15
  qtable[4,]<-c(0.300, 0.342, 0.425)  #size=20
  qtable[5,]<-c(0.277, 0.317, 0.393)  #size=25
  qtable[6,]<-c(0.260, 0.298, 0.372)  #size=30

  nn<-length(x)

  qstat.low<-(x[rank(x)==2]-x[rank(x)==1])/(max(x)-min(x))
    qstat.high<-(x[rank(x)==nn]-x[rank(x)==(nn-1)])/(max(x)-min(x))

  qstat<-ifelse(qstat.low>qstat.high, qstat.low, qstat.high)
  ind<-ifelse(qstat.low>qstat.high,(1:length(x))[rank(x)==1], (1:length(x))[rank(x)==nn])

print(x)
  cat("qstat.low=", qstat.low, "\n")
    cat("qstat.high=", qstat.high, "\n")
    cat("qstat=", qstat, "\n")
    cat("ind=", ind, "\n")

  if (prob==0.90) k<-1
  if (prob==0.95) k<-2
  if (prob==0.99) k<-3

 ind0<-rep(FALSE, size)

if (size>=5) {
  if (qstat>qtable[min(floor(size/5),6),k])
   ind0[ind]<-TRUE
}
    return(ind0)
}



out.test<-function(x, size, prob, method="ttest") {

if (method=="ttest") {
  x0<-abs(x-mean(x))/(var(x)^0.5)
  prob0<-1-pt(x0, size-1)
  ind<-(1:length(x0))[prob0<(1-prob)/2]
  return(list(outliers=any(prob0<(1-prob)/2), ind=ind))
}
if (method=="qtest")
  return(qtest(x,size,prob))
}

out.test<-function(x, size, prob, method="ttest") {

if (method=="ttest") {
  x0<-abs(x-mean(x))/(var(x)^0.5)
  prob0<-1-pt(x0, size-1)
  ind<-(1:length(x0))[prob0<(1-prob)/2]
  return(list(outliers=any(prob0<(1-prob)/2), ind=ind))
}
if (method=="qtest")
  return(qtest(x,size,prob))
}
