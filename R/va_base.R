## this program implements the method discussed in King and Lu(2008)
## details see Appendix of King & Lu (2008)

# ???what do I need in Design?

require(Design, quietly=TRUE)

# We use quadprog:solve.QP in quad.constrain.R
# for constrained linear optimization
require(quadprog, quietly=TRUE)


## va() is an interface function that calls cod.est.base()
## the purpose of va() is to make inputs more user friendly and check data
## It doesn't work well if the data is not cleaned.
## In a Cleaned data: by column there are should be some variability,
## if all responses=1 or 0,then delete


##main function
## for parameters, see ?va

##??? do I need checkrow.names?
va<-function(formula, data=list(hospital=NA, community=NA),nsymp=16,
             n.subset=300, method="quadOpt", fix=NA, bound=NA,
             prob.wt=1,  boot.se=FALSE,
             nboot=300, printit=TRUE,
             print.reg.size=TRUE, checkrow.names=FALSE,
             predict.S=FALSE){

  # if number of bootstrap=1, will not computer bootstrapped
  # based standard errors
  if (nboot==1) boot.se<-FALSE

  ## check if quadprog is loaded
  if (method=="quadOpt") {
    tt<-requireNamespace("quadprog")
    if (!tt)
      stop("\n Require Quadractic Programming package 'quadprog'\n",
           call.=FALSE)
  }

  # users can apply prior knowledge to the estimation
  # prior can be two types: fixed value, or bounds
  # see manual (fix and bound) for details.
  prior <- list()

  # fFix and fBdd are two internal functions that reads
  # information from inputs "fix" and "bound" and assign
  # the values to a list object "prior"
  # definitions of fFix and fBdd are at the end of this file

  if (!is.na(fix)) {
    prior$fix <- fFix(fix)
  }
  if (!is.na(bound)) {

    prior$bdd <- fBdd(bound)
  }

  # offers users an opportunity if priors
  # are read correctly
  # An interactive interface would be better, for example
  # "You have input priors :....
  #  to proceed, hit enter
  #  to modify the priors, type 0 "

  print(prior)


  # all.vars() and dot() are internal functions
  # that extract colnames from the formula
  # definitions of all.vars and dot are at the end of this file
  rhsVars<-all.vars(formula[[3]])
  lhsVars<-all.vars(formula[[2]])
  lhsVars<-dot(lhsVars,data[[1]])


  # rhsVars contains colnames corresponding to
  # the cause of death
  # We allow the cause of death variable to have different
  # colname in the hospital and community samples
  # or missing in the community(population) sample
  if(length(rhsVars)==2){
    hospVars<-c(rhsVars[[1]],lhsVars)
    commVars<-c(rhsVars[[2]],lhsVars)
  }
  else
    hospVars<-commVars<-c(rhsVars,lhsVars)


  ## first column cause of death
  ## second to last symptoms
  hospD <- data[[1]][,hospVars] #key data matrix

  ## Sometimes cause of death variable (rhs[[1]]) is not present
  ## in the community dataset. In this case, add one variable with
  ## the same name and arbitrary numeric values (read documentation)

  comm.cod <- TRUE
  commD <- data[[2]][,commVars]
  if ((rhsVars[[1]] %in% colnames(commD))==FALSE){
    warning("cause of death not in population sample")
    comm.cod <- FALSE
    tmp <- rep(1,nrow(commD))
    commD<-cbind(tmp,commD)    #key data matrix
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


  # obtain a complete list of the causes of death
  # note if a cod is not present in the hospital sample,
  # the estimate will be 0 in the communit sample.
  cod.names<-union(names(table(hospD[,1])), names(table(commD[,1])))

  index.all<-cod.names
  index.all<- index.all[order(as.numeric(index.all))]
  # number of cods
  nd <- length(index.all)

  # number of symptoms
  ns <- ncol(hospD)-1

  # if user does not supply prob.wt
  # prob.wt=1  binomial weights are assinged
  # prob.wt=0 no weights assigned
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


  ######################
  ## Key estimation step
  ## cod.est.base estimate the COD in community sample
  ## based on  'n.subset' numbre of subsets, each of which has
  ## randomly (with prob.wt) selected  nysmp number of symptoms

    res0<- cod.est.base(hospital=hospD,community=commD,
                         n.subset=n.subset, nsymp=nsymp,
                         method=method,prior=prior,
                         prob.wt=prob.wt,printit=printit,
                         print.reg.size=print.reg.size,
                         checkrow.names=checkrow.names)


  ## calculate s.e. only when boot.se and nboot>1
  if (!boot.se) nboot<-1

  if (nboot==1) boot.se<-FALSE


  # comp is a nboot by nd by 2 array saving the bootstrap results
  # if in community sample, true cod is know, it will be saved under "true"
  # "emp" saves the estimated cod based on each bootstrap sample.
  if (boot.se) {
    comp<-array(0, dim=c(nboot, nd,2),
              dimnames=list(1:nboot, index.all, c("true", "emp")))

  # s.fit will save the predicted P(S) for each bootstrapped
   # community sample
    s.fit <- matrix(0, nboot, ncol(hospD)-1)


  # start bootstrapping....
    for (t in 1:nboot){
      p1<-dim(commD)[1]
      p2<-dim(hospD)[1]

      boot.good<-FALSE

      while (!boot.good && boot.se) {

        # randomly sample row indices with replacement
        index.boot1<-sample(1:p1, p1, replace=TRUE)
        index.boot2<-sample(1:p2, p2, replace=TRUE)

        # bootstrapped samples
        commD0<-commD[index.boot1,]
        hospD0<-hospD[index.boot2,]

        # if the bootstrapped sample of either samples
        # has any sympotms with variation, return to
        # draw another bootstrapped samples.
        boot.good<-!(any(colMeans(commD0[,-1])==1) ||
                     any(colMeans(commD0[,-1])==0) ||
                     any(colMeans(hospD0[,-1])==1) ||
                     any(colMeans(hospD0[,-1])==0) )
      }


      # estimate COD baed on bootstrapped samples
      res.boot<- cod.est.base(hospital=hospD0,community=commD0,
                         n.subset=n.subset, nsymp=nsymp,
                         method=method,prior=prior,
                         prob.wt=prob.wt,printit=printit,
                         print.reg.size=print.reg.size,
                         checkrow.names=checkrow.names)

      cod.emp<-res.boot$est.CSMF
      cod.true<-res.boot$true.CSMF

      # store estimates and "true"(could be junk) in array "comp"
      comp[t,index.all %in% row.names(cod.true), 1]<-cod.true
      comp[t,index.all %in% names(cod.emp), 2]<-cod.emp

      # predict s.fit, or, P(S) for each boortrapped community sample
      # this can be used to calculate s.e of P(S)

      if (predict.S==TRUE) {
        nd <- length(index.all)
        ns <- ncol(hosp)-1
        sd.mat<-matrix(0,ns,nd)
        for (k in 1:nd)
          {
            sd.mat[,k] <- colMeans(hospD0[hospD0[,1]==k,2:(ns+1)])
          }
        sd.mat[is.na(sd.mat)]<-0
        # s.fit is clculated by P^h(S|D)*P(D)
        # p^h(S|D) is estimated based on a bootstrapped sample of hospital
        # p(D) is the COD estimates based on one sets of bootrapped samples
        s.fit[t,]<-sd.mat%*%res$est.CSMF

      }


    cat(paste(t, "th bootstrap sampling done.\n\n"))
    }
  }  # end of bootstrap



    ### outputs

    out<-list()

    out$est.CSMF<-rep(0,nd)
    out$true.CSMF<-rep(0,nd)
   # out$est.CSMF saves the COD estimates
   # the names of the elements of vector are union of hospital and
   # community COD.
  out$est.CSMF[index.all %in% names(res0$est.CSMF)]<-res0$est.CSMF
    names(out$est.CSMF)<-index.all


  # out$true.CSMF saves the true COD in the community if available
    if (comm.cod) {
      out$true.CSMF[index.all %in% names(res0$true.CSMF)]<-res0$true.CSMF
      names(out$true.CSMF)<-index.all
    }
    else
      out$true.CSMF<-NA

   # estimated COD based on each subset of syptoms
   # can be used to check quality of estimates
   # an unstable one could have many 0s in any set of COD estimate
    out$subsets.est<-res0$est.res

  #? what is index.match?
    if (n.subset==1)
      out$index.match<-res0$index.match


   if (boot.se) {
        # estimated s.e. of out$est.CSMF
        out$CSMF.se<-apply(comp[,,2], 2, var)^0.5
        # bootstrap mean CSMF of the community sample if
        # truth is known
        out$true.CSMF.bootmean<-colMeans(comp[,,1])
        # estimated s.e. of true COD
        # can be used as a benchmark s.e.in the COD distribution
        out$true.bootse<-apply(comp[,,1], 2, var)^0.5
        # estimation result matrix nboot*nd*2
        out$res.boot<-comp

        # output bootstrap results for fitted P(S)
        if (predict.S==TRUE) {
          out$s.boot.fit<-s.fit
        }
      }
    class(out)<-"VA"
    return(out)
  }

## Key function  implementing quadratic optimization
## or constrain least square based on subsets of symptoms
##hospital is the hospital based data  ???
##community is the community based data  ???
##format--MRcause, Symptom 1,2,...

cod.est.base<-function(hospital,community, n.subset=300, nsymp=16,
                       prior=list(), prob.wt=NULL, printit=TRUE,
                       print.reg.size=TRUE, method=NULL,
                       checkrow.names=FALSE, h.cod=NULL){

        # default method is "quadOpt"
        if (is.null(method)) method<-"quadOpt"

        # q equals to number of available symptoms +1
        q<-dim(community)[2]

        # if the size of the subset is larger than available
        # number of symptoms, produce an error message
        # user should pick a smaller nsymp
        if (nsymp > (q-1))
          stop("nsymp is bigger than actual no. of symptoms")

        # the number of subsets need to be smaller than
        # the number of different subsets of nsymp symptoms
        # that can be drawn from the total (q-1) sysmptoms
        # otherwise reutrn a warning
        if (n.subset>choose((q-1), nsymp))
          warning("n.subset is bigger than choose(q, nsymp)")

        # total number of COD in the hospital sample
        nd<-length(table(hospital[,1]))

        # d.res is an internal matrix that helps to monitor the code
        # remove d.res this eventually
        d.res<-matrix(NA, n.subset, nd)
        colnames(d.res)<-names(table(hospital[,1]))

        # start subsetting...
        for (subset in 1:n.subset){
                sample.ok<-FALSE
                counter<-1
                while (!sample.ok) {
                        sample.ok<-TRUE
                        # subset symptoms
                        index<-sample(2:q, nsymp, replace=FALSE, prob=prob.wt)


                        # h.probs contains P(S) and P(S|D) in hospital smaple
                          h.probs<-DISTgen(hospital[,c(1,index)], whichdata="hospital")
                        # c.probs contains P(S) and P(S|D) in community smaple
                        c.probs<-DISTgen(community[, c(1,index)],whichdata="community")

                        # Y=X beta
                        # Y=P(S), X=P(S|D), beta=COD to be solved

                          # P(S)
                          y<-c.probs$tsp
                          # P^h(S|D)
                          x<-h.probs$sd.mat

                        #match the S components in both quantities.
                        # so Y and X have same rows

                          index.match<-intersect(row.names(y), row.names(x))

                      # y.use and x.use are subsets of P(S) and P^H(S|D)
                      # that have same sysmptom profiles
                      # note S is S_{sub} only

                         y.use<-y[row.names(y) %in% index.match]
                          x.ind<-row.names(x) %in% index.match

                        # P(S|D) tend to be sparser
                          x.use<-NULL
                          if (!is.null(x.ind))
                           x.use<-x[x.ind,]


                        # not use
                          ratio.y<-length(y.use)/length(y)
                        # not use
                        ratio.x<-nrow(x)/nrow(x.use)

                        ## y.use and x.use need to be non empty
                        ## nrow in x.use needs to be larger than ncol of x.use

                          if (is.null(y.use) || is.null(x.use)
                              || is.null(dim(x.use))
                              || (!is.null(dim(x.use)) &&
                                  (dim(x.use)[1]<dim(x.use)[2]))
                              ) {

                            # if above condition is not met,
                            # go back to take a different smaple

                            sample.ok<-FALSE
                            counter<-counter+1



        # one main reason cause !sample.ok is due to sparsity
        # if sample.ok appears to be a problem more than twice
        # output a warning
           # if sample.ok appears to be a problem more than ten times
        # stop and suggest user to pick smaller nysmp (so P(S|D) can be
                            #less sparse


                            if (counter==2)
                              warning("sample could be too low
                                  for nsym\n")
                                  else
                                    if (counter==10)
                                      stop("please use a smaller nsymp\n", call.=FALSE)

                          }
                        if (method=="quadOpt") {
                          # estimate COD using quad.constrain
                          # if no solution, coef return "0"
                          coef<-try(quad.constrain(y.use, x.use, prior=prior), TRUE)
                          # if error occurs, go back to sample another subset
                          if (length(coef)==1)
                            sample.ok<-FALSE
                        }

                      }
                # user can monitor the sparsity of the X matrix, P(S|D) or x.use
                if (print.reg.size) {
                  cat("the dimension of X matrix in constraint
                optimization step\n")
                  print(dim(x.use))
                }


                # estimate COD using ls.constrain
                # note this method is not sensitive on data input
                # it will estimate CODs regardless the sparsity of x.use
                if (method=="constrainLS")
                  coef<-ls.constrain(y.use, x.use, sum1=TRUE, prior=prior)

                # save estimates for each subset
                d.res[subset, 1:nd]<-coef

                # user can monitor the progress of hte code
                if ((printit) & (subset==round(subset/n.subset,digits=1)*n.subset))
                  cat(paste(round(subset/n.subset, digits=1)*100, "percent done. \n"))

              }

        # cod is the true COD in the community, if available
        # otherwise junk
        cod<-table(community[,1])
        cod<-cod/sum(cod)
        c.probs$cod<-cod

        # index.match is an internal object
        # I used it to test the code
        if (!checkrow.names) index.match<-NULL


        if (n.subset>1){
          # make sure Sum=1 constrain is respected
          good<-rowSums(d.res[,1:nd])==1
          est<-colMeans(as.matrix(d.res[good,1:nd]), na.rm=TRUE)
        }
        else if (n.subset==1)
          est<-coef
        res<-list(est.CSMF=est, true.CSMF=c.probs$cod,
                  est.res=d.res, index.match=index.match)

        return(res)

      }



dot<-function(lhsVars,data){
        dotpos<-match("...",lhsVars)
        if(is.na(dotpos))
          return(lhsVars)
        tmp<-colnames(data)
        if(dotpos==1){
                if(length(lhsVars)==1)
                  return(tmp)
                out<-tmp[1:which(tmp==lhsVars[[dotpos+1]])-1]
                lhsVars<-lhsVars[-(1:dotpos)]
                dotpos<-match("...",lhsVars)
        }else
        out<-c()
        while (!is.na(dotpos)){
                bef.miss.index<-which(tmp==lhsVars[[dotpos-1]])
                if(dotpos==length(lhsVars))
                  after.miss.index<-length(tmp)+1
                else
                  after.miss.index<-which(tmp==lhsVars[[dotpos+1]])
                out<-c(out,lhsVars[1:(dotpos-1)],tmp[(bef.miss.index+1):(after.miss.index-1)])
                lhsVars<-lhsVars[-(1:dotpos)]
                dotpos<-match("...",lhsVars)
        }
        out<-c(out,lhsVars)
        return(out)
      }




##pack fix and bound input into list like prior

fFix <- function(fix) {

     getNamesF<-function(x)
          x[[1]]
     getValuesF <- function(x)
          x[[2]]

    parseF <- function(x)
          sapply(x, strsplit,split="=")

        pfixed <- sapply(fix,parseF)
        fix <- matrix(as.numeric(sapply(pfixed,getValuesF)),
                        dimnames=list(sapply(pfixed,getNamesF),"value"))
     fix
   }

fBdd <- function(bound) {

        getNamesB <- function(x)
          x[[2]]
        getValuesB <- function(x)
          c(x[[1]],x[[3]])

        parseB <- function(x)
          sapply(x,strsplit,split ="<")

        pbound <- sapply(bound, parseB)
        bdd <- matrix(as.numeric(unlist(sapply(pbound, getValuesB,
simplify=FALSE))),
                      nrow=length(pbound), dimnames=list(sapply(pbound,
getNamesB),c("lbd","ubd")))

        bdd
}
