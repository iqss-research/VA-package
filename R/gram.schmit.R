gram.schmidt<-function(mat, orthnorm=NULL)
  {
    if (det(mat)==0) stop("mat is not full rank")

    omat<-mat
    p<-dim(omat)[1]

    for (i in 2:p)
      {
        temp <- rep(0,p)
        for (j in 1:(i-1))
          temp <- temp+t(omat[j,])%*%
      omat[i,]/(t(omat[j,])%*%omat[j,])*omat[j,]
        omat[i,] <- omat[i,]-temp
      }

    ##orthornomalize certain row of omat
    ## if orthnorm=NULL, no normalization
    ## if orthnorm=vector, index whic rows to be normalized
    
    if (!is.null(orthnorm)) {
        oind<-orthnorm
        for ( i in oind)
        omat[i,] <- omat[i,]/(t(omat[i,]%*%omat[i,]))^0.5
    }
    return(omat)

  }
