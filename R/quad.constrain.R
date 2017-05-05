###Constrained Least Squares to estimate \sum(beta)=1, beta>0

quad.constrain<-function(Y, X,  prior=list())
  {
    p<-dim(X)[1]
    q<-dim(X)[2]

    ###prior is a list that takes two kinds of contrains
    ## not tested yet one is fix value, the other is interval constrain
    ind <- matrix(0,q,2)
    rownames(ind) <- colnames(X)
    # if ind[,1]=0, cod fixed values,
    # =1, cod to be stimated 
    ind[,1] <- 1

    
    const<-1
    Y.new<-Y
    if (!is.null(prior$fix)) {
        for (i in 1:dim(prior$fix)[1]) {
          id<-(1:q)[rownames(prior$fix)[i]==rownames(ind)]
          value<-prior$fix[i]
           ind[id,1]<-0
           ind[id,2]<-value

           const<-const-value

           Y.new<-Y.new-X[,id]*value
           if (const<0) stop("invalid inputs of known cods.\n")    
        }
      }

    

    lbd <- rep(0,q)
    ubd <- rep(1,q)
    if (!is.null(prior$bdd)) {
      for (i in 1:dim(prior$bdd)[1]) {
        id <- (1:q)[rownames(prior$bdd)[i]==rownames(ind)]
        lbd[id] <- prior$bdd[i,1]
        ubd[id] <- prior$bdd[i,2]
      }
    }
    
    X0 <- X[,ind[,1]==1]
 
    Dmat<-t(X0)%*%X0
 
    dvec<-t(Y.new)%*%X0
    q0<-ncol(X0)
    Amat<-matrix(0, q0,q0*2+1)
    Amat[,1]<-rep(1,q0)
    Amat[,2:(q0+1)]<-diag(1,q0)
    Amat[,(q0+2):(2*q0+1)]<-diag(-1,q0)
    
#   Amat<-t(rbind(rep(1,q0), diag(q0), -diag(q0)))
   bvec<-c(const, lbd[ind[,1]==1], -ubd[ind[,1]==1])

    # minimizing (Y-x\beta)^2 subject to sum(beta)=const

    
    
   res<-solve.QP(Dmat,dvec, Amat, bvec, meq=1)$solution

    ind[(ind[,1]==1),2] <- res
   
    return(ind[,2])
    
  }


