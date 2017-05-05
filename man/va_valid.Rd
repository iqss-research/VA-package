\name{va.validate}
\alias{va.validate}

\title{va.valid: a method selecting potentially biased symptoms}
\description{
  This function selects the symptoms that are potentially reported with
  bias by comparing the fitted marginal distribution of the symptoms
  with the observed symptom distribution in the community sample. To
  control for the false discovery rate, Bonferroni adjustment to the
  significance values can be sequentially applied to the models.
}
  

\usage{ va.validate(formula,  data=list(hospital=NA, community=NA),
               nsymp=10,n.subset=300, nboot=1, boot.se=FALSE, method="quadOpt",
                      fix=NA, bound=NA, prob.wt=1, 
                      printit=TRUE,print.reg.size=TRUE,
                      clean.method="ttest",
                      min.Symp=10, confidence=0.95, FDR=TRUE)
}

\arguments{
  \item{formula}{A formula object. The left side of the formula is
    the collection of symptoms. The right side is the cause of
    death. For example, if there are 5 symptoms, named
    \code{fever},\code{coughing},\code{chestpain},\code{dizziness},
    \code{shortbreath}, and the cause of death variable is \code{death},
    then the formula can be written as:

    \code{formula=cbind(fever, coughing, chestpain, dizziness,
      shortbreath)~death}

    or for short as:
    \code{formula=cbind(fever, ... ,shortbreath)~death}

    Note that the short way of writing formula requires the symptoms variables
    are located in a consecutive block in the data starting from
    \code{fever} and ending with \code{shortbreath}.
   Note that the current version requires the varible on the right hand
   side of the formula, \texttt{death} in this example, to be present in
   the \texttt{community} sample. If it is unknown in the
   \texttt{community} sample, the user needs to create such variable
   with arbitrary numerical values. }
  \item{data}{a list of two datasets. The first is the hospital
    data, which contains the known cause of death for each
    individual, and a collection of symptoms from verbal autopsy studies.
    The second is the community data where typically only the symptoms are
    available. The known
    cause of death can be available outside hospital
    if it is a validation study, but it will not be used during
    estimation. Variable names must be exactly the same in two data
    sets.}
 \item{nsymp}{a positive integer, specifing the size of subsets of
    symptoms drawn from the total set for estimating cause specific
    mortality fractions at each iteration. The optimal number of
    \code{nsymp} can be found calling \code{va_gcv}, which use
    general cross-validation method to find the optimal size of
    subset that minimize the prediction
    errors based on the training data(typically, hospital data).
    For more details, refer to King and Lu (2008).}
  \item{n.subset}{A positive integer specifing the total number of
     subsets and thus estimations of all symptoms.
     The default is \code{300}.}
  \item{nboot}{a positive integer. If \code{boot.se=TRUE}, it
    specifies the number of bootstrapping samples taken to estimate the
    standard errors of CSMF. The default is \code{1}.}
   \item{boot.se}{a Logical value. If \code{TRUE}, bootstrap
     standard errors of the CSMF are estimated.  This typically takes a lot
     of computing time. It is highly suggested to set
     \code{boot.se=FALSE} in \code{va_gcv}. Default=\code{FALSE}.}
      \item{method}{ A string specifying the computational procedure
     used to estimate the cause specific mortality fractions. When 
     \code{method=''quadOpt''}, CSMF is estimated via constrained 
     quadratic programming. A subroutine (\code{Solve.QP}) from the 
     \code{quadprog} package is called to perform the constrained 
     quadratic optimization task. When \code{method=``constrainLS''},
     CSMF is estimated via constrained least squares. The default method
     is \code{quadprog} as it is faster and more stable.}
   \item{fix}{ A vector of strings that specifies whether a subset of
     the cause specific mortality fractions are set to predetermined
     values (based on, e.g.,the information obtained from other
     sources). Suppose we would like to prefix ''d1'' to be 5\%, ''d2''
     to be 15\%, then \code{fix=c("d1=0.05", "d2=0.15")}. The default
     is \code{NA}, no such constrain is imposed. }
   \item{bound}{ A vector of strings that specifies lower and upper
     bounds of a subset of the cause specific mortality fractions
     (based on, e.g.,the information obtained from other
     sources). Suppose we would like ''d3'' to be estimated between 5\% and
     10\%, "d4" to be between 1\% and 2\%, then
     \code{bound=c("0.05<d3<0.1", "0.01<d4<0.02")}. The default
     is \code{NA}, no such constrain is imposed. }
   \item{prob.wt}{A positive integer or a vector of weights that determines how
     likely a symptom is of being selected for a subset. When
     \code{prob.wt} is a user input vector, it needs to be a vector 
     of probabilities and sum up to 1. The length of
     \code{prob.wt} needs to be equal to the total number of symptoms.
     When \code{prob.wt=1}, binomial weights which are proportion to 
     the inverse of variances of the each reported binary symptom variable.   
     When \code{prob.wt=0}, all symptoms will be equally selected. 
     The default is \code{1}.}    
  \item{printit}{Logical value. If \code{TRUE}, the progress of the
    estimation procedure will be printed on the screen.}
  \item{print.reg.size}{Logical value. If \code{TRUE}, the size of the
  regression matrix is printed at each step of subsampling. It provides
  helpful information for user to choose the number of symptoms to
  subsample. It is recommended to print the size of the regression
  matrix for different values of \texttt{nsymp} with a small size of
  \texttt{n.subset}.}
\item{clean.method}{A string specifying which test to use to detect
    poorly fit symptoms. The default is \code{''ttest''}. The other
    option is \code{''ztest''}. For details, see King and Lu (2008b).}
  \item{min.Symp}{An integer value. When the number of availability
    symptoms is less than \code{min.Symp}, the automatic procedure for
    selecting next poorly fit symptom will stop. The default is 10.
    But if \code{min.Symp} is less than \code{nsymp}, \code{va.validate}
    coerces \code{min.Symp} to be \code{nysmp}.}
  \item{confidence}{A number between 0 and 1. It specifies the
    confidence level (or the significance level \code{(1-confidence)/2})
    at which user decided whether to remove a symptom out of the
    estimation. The default is 0.9
  }
  \item{FDR}{Logical value. If \code{TRUE}, a Bonferroni adjustment for
    multiple testing is applied sequentially to a collection of nested
    models as more symptoms being removed. For details, see King and Lu (2008b).} 
}


\details{
  For details, please refer to ''Designing Verbal Autopsy Analyses: A
  Report to the World Health Organization'' (King and Lu, 2008b) and
  \url{http:\\gking.harvard.edu\va} 
}

\value{
  \code{va.validate} outputs the following objects. \code{cod.list}
  returns a list of cause of death estimations based on a set of nested
  models. \code{Ps.list} returns a list of fitted marginal symptom
  distributions based on the same set of nested
  models. \code{delete.list} returns a list of collections of removed
  symptoms based on the nested models. \code{FDR.delete.list} returns a
  list of symptoms that should be removed based on the Bonferroni
  adjustment and attained a global confidence level of
  \code{confidence}. When \code{boot.se} is \code{TRUE},
  \code{va.validate} also returns a set of objects that summarizes the
  variance of the predicted marginal symptom distribution
  (\code{Ps.var}), the variance of the residuals (\code{e.var}), the
  results based on all the bootstrapped samples (\code{cod.list.boot}
  and \cod{Ps.list.boot}).
}
 

\references{
  King, Gary and Ying Lu. (2008) \dQuote{Verbal Autopsy Methods with
    Multiple Causes of Death}, 14(1), Statistical Science. Also
  available at \url{http:gking.harvard.edu/va}
  King, Gary and Ying Lu. (2008b) \dQuote{Designing Verbal Autopsy
  Analyses: A Report to WHO}. 

}

\keyword{verbal autopsy, CSMF}
  
