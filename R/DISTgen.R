## make p(S), P(D) and p(S|D) from the hospital data
## with option of boostraping symptoms
## wmat: Data takes the format
## column 1: causes of disease
## column 2--q: symptoms


DISTgen<-function(wmat, whichdata) {
  
        p<-dim(wmat)[1]
        q<-dim(wmat)[2]

        tsp<-sd.mat<-NULL

        # if fewer than 42 symptoms, treat any row of wmat
        # as binary sequence, 010001....1100
        # and represent it as decimals.
       if (q<43) {
          sp <- apply(wmat[,2:q],1,bintodec)
        }else{
          ##otherswise, if more than 42 symptoms,
          ## treat "010010001...0011" as string
          sp <- apply(wmat[,2:q],1,paste,sep="",collapse="")
        }
        # calculate P(S) by tabulation in the community
        if (whichdata=="community") {
           tsp <- table(sp)
           tsp<-tsp/sum(tsp)
         }

        # caluclate P(S|D) by tabulation in the hospital
        if (whichdata=="hospital") {
          sd.mat <- table(sp, by=wmat[,1])
          for (j in 1:ncol(sd.mat))
            sd.mat[,j]<-sd.mat[,j]/sum(sd.mat[,j])
        }

        # calculate cod in either sample by tabulating COD.        
        cod <- table(wmat[,1])        

        cod<-cod/sum(cod)
 
        return(prob.dist=list(tsp=tsp, sd.mat=sd.mat, cod=cod))
}

 bintodec<-function(x){
 sum(x * 2^ ((length(x)-1):0))
 }
