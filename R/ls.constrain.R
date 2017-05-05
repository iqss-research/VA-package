###Constrained Least Squares to estimate \sum(beta)=1, beta>0
## See appendix

ls.constrain<-function(Y, X, sum1=TRUE, prior=list(), greater.zero=TRUE)
  {
    # Number of obs (rows of P(S|D) and P(S))
    p<-dim(X)[1]
    # Number of CODs to be estimated
    # length of beta
    q<-dim(X)[2]

    ###prior is a list that takes two kinds of contrains
    ##one is fix value, the other is interval constrain

    # ind is a matrix that keeps track of
    # the parameter estimation of beta (i.e. P(D), cod)
    ind<-matrix(0, q, 4)
    # each row of ind corresponds to a cause of death
    # i.e. coefficient id
    rownames(ind)<-colnames(X)

    # firt column of ind indicate if a beta is fixed
    # through estimation
    ind[,1]<-1     #fixed or not, 0 being fixed
    const<-1
    Y.new<-Y

    # attach prior knowledege to ind matrix
    if (!is.null(prior$fix)) {
        for (i in 1:dim(prior$fix)[1]) {
           # check which betas are specified to be fixed
           id<-(1:q)[rownames(prior$fix)[i]==rownames(ind)]
           value<-prior$fix[i,1]
           # fix those beta
           ind[id,1]<-0
           # to the values specified in prior$fix
           ind[id,2]<-value
           # ??? what is 3rd column for
           ind[id,3]<-0
           # once a beta is fixed to "value", the sum_to_1 constrain
           # is reduced by "value"
           const<-const-value
           # Y=b1*X1+b2*X2+b3*X3
           # Y.new=Y-b10*X1, suppose b1 is fixed to b10
           Y.new<-Y.new-X[,id]*value
           if (const<0) stop("invalid inputs of known cods.\n")    
    }
   }

    ##check if there is enough information, 
    ## if too many coefficients are NA, then skip

    # right now this step is done without incorporating
    # any prior information: inconsistency?

    # LS solution of beta
    coef0<-lm(Y~X-1)[[1]]
   pass<-FALSE
    # if a simple LS without constraints can't
    # estimate more than 2/3 of hte coefficient
    # then the data is of poor quality

    ind[is.na(coef0),1]<-9
    



    # number of COD that needs to be estimated (not fixed)
    nd<-length(ind[ind[,1]==1,1])


    if (nd<floor(q/3))
      {
        pass<-TRUE
        ind[,2]<-NA
	  ind[,4]<-NA
      }

    while (!pass & (nd>=floor(q/3)) & nd>=2) {
       pass<-TRUE
       ##constraint C'beta=1, C={1/sum(beta)}={1/d}

       # G is an orthonormal matrix
     
       G<-diag(1,nd)
        G[1,]<-rep(1/const,nd)
        G<-as.matrix(gram.schmidt(G, orthnorm=2:nd))
        C<-G[1,]
        A<-G[2:nd,]
         ##decompose X=ZA+WC=(W,Z)%*%G
         ##(W,Z)=Xsolve(G), G=(C',A')'
	   mtemp<-X[,ind[,1]==1]%*%solve(G)
         Z<-mtemp[,2:nd]
         W<-mtemp[,1]

         #Y=X*beta=Z*A*beta+W1'*beta=Z*A*beta+W
         #Y-W=Z*gamma, gamma=A*beta
         Y.use <- Y.new-W

         ##G*beta=(c*beta, A*beta)=(1, gamma)
         ##beta=solve(G)*gamma
         ##cov(beta)=solve(G)%*%[vcov(gamma)+ZERO(1,1)]%*%t(solve(G))
         ##ZERO[1,1], 1st row and 1st column are 0
         reg <- lm(Y.use~Z-1)
	   gamma<-reg[[1]]
         gamma.cov<-vcov(reg)
 	  # invZ<-ginv(t(Z)%*%Z)
        # gamma<-invZ%*%t(Z)%*%Y.use
     	  # sigma2<-as.numeric((t(Y.use)%*%Y.use- t(Y.use)%*%Z%*%invZ%*%t(Z)%*%Y.use)/(length(Y)-nd))
        # gamma.cov<-invZ*sigma2
         beta <- solve(G)%*%c(1,gamma)
         gamma0.cov<-matrix(0,nd,nd)
         gamma0.cov[2:nd,2:nd]<-gamma.cov
         beta.cov<-solve(G)%*%gamma0.cov%*%t(solve(G))
         tval<-abs(beta/sqrt(diag(beta.cov)))
       # beta 
       ind[(ind[,1]==1),2]<-beta
       # tvalue of beta, test whether beta<0
       ind[(ind[,1]==1),3]<-tval
       # standrad error of beta
         ind[(ind[,1]==1),4]<-sqrt(diag(beta.cov))
      
	   if (greater.zero) {
          if (any(beta<0))
            {
              # set the negative beta with largest negative t-value to be 0
              ind[(ind[,2]<0) & (ind[,3]==max(ind[(ind[,2]<0),3])),1:3]<-0
              nd<-length(ind[ind[,1]==1,1])
        # if hte number of nonzero betas are less than 1/3 of nd, stop
              if (nd<floor(q/3)) ind[,3]<-NA
              pass<-FALSE
            }
          ###build the interval constraint here.
	  }

   }
    sse<-sum((Y-as.vector(X%*%ind[,2]))^2, na.rm=TRUE)
    return(ind[,2])
    # need to return a message if the optimization is not completely done
    # when quit at nd<q/3
  }



###not used..
ls.constrain.old<-function(Y, X, sum1=TRUE, greater.zero=TRUE, prior=list())
  {
    p<-dim(X)[1]
    q<-dim(X)[2]

    ind<-matrix(0, q, 4)
    ind[,1]<-1:q
    ind[,2]<-1
    coef0<-lm(Y~X-1)[[1]]
    ind[is.na(coef0),2]<-9
    pass<-FALSE


    nd<-length(ind[ind[,2]==1,1])
    if (nd<floor(q/3))
      {
        pass<-TRUE
        ind[,3]<-NA
      }

    while (!pass & (nd>=floor(q/3))) {
       pass<-TRUE
      if (sum1) {
        omat<-diag(1,nd)
        omat[1,]<-rep(1,nd)
        omat<-as.matrix(gram.schmidt(omat, orth=2:nd))
      
        omat2<-matrix(0, nd-1, nd-1)
        for ( i in 1:(nd-1))
        omat2[,i] <- omat[2:nd,i]-omat[2:nd,nd]
      }
      
      if (!sum1) {
        coef<-lm(Y~X[,ind[(ind[,2]==1),1]]-1)
      tval<-abs(summary(coef)$coefficients[,3])
      ind[(ind[,2]==1),3]<-coef[[1]]
      ind[(ind[,2]==1),4]<-tval
      
      if (any(coef[[1]]<0))
        {
          ind[(ind[,3]<0) & (ind[,4]==max(ind[(ind[,3]<0),4])),2:4]<-0
          pass<-FALSE
        }
    }
      if (sum1)
        {
          if (nd>2)
            {
              X.use <- X[,ind[(ind[,2]==1),1]]%*%t(omat[2:nd,])
            }
          if (nd<=2)
            {
              X.use <- X[,ind[(ind[,2]==1),1]]%*%omat[2:nd,]
            }

          W <- (X[,ind[(ind[,2]==1),1]]-X.use%*%omat[2:nd,])[,1]
          Y.use <- Y-W

          coef2 <- lm(Y.use~X.use-1)
          coef <- solve(omat2)%*%(coef2[[1]]-omat[2:nd,nd])
          cov.coef<-solve(omat2)%*%vcov(coef2)%*%t(solve(omat2))
          tval<-abs(coef/sqrt(diag(cov.coef)))
          coef<-c(coef, 1-sum(coef))
          tval<-c(tval,0)
          ind[(ind[,2]==1),3]<-coef
          ind[(ind[,2]==1),4]<-tval

          if (any(coef<0))
            {
              ind[(ind[,3]<0) & (ind[,4]==max(ind[(ind[,3]<0),4])),2:4]<-0
              nd<-length(ind[ind[,2]==1,1])
              if (nd<floor(q/3)) ind[,3]<-NA
              pass<-FALSE
            }
          
        }
    }
    return(ind[,3])
    
  }


