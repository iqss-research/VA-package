\name{va}
\alias{va}

\title{VA:Software for Analyzing Verbal Autopsy Data}
\description{
  Estimates cause-specific mortality rates in a population where a
  set of dichotomous symptoms are available, using the relationship 
  between symptoms and a  multicategory cause-of-death
  variable collected from a nearby medical facility.  Estimation is
  nonparametric.}

\usage{va(formula, data=list(hospital=NA,community=NA), nsymp=16, n.subset=300,
method=``quadOpt'', fix=NA, bound=NA, prob.wt=1, boot.se=FALSE,
nboot=300, printit=TRUE, print.reg.size=FALSE, predict.S=FALSE)
}

\arguments{
  \item{formula}{a formula object. The left side of the formula is
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
    mortality fractions at each iteration. \code{nsymp} can be found
    calling \code{va_gcv}, which use general cross-validation method
    to find the optimal size of subset that minimize the prediction
    errors based on the training data(typically, hospital data).
    For more details, refer to King and Lu (2006).
    For practical purpose, we give the following 
  recommendataions: for total number of causes of death \code{D<=10}, 
  use 7-12 symptoms; for \code{D>10}, use 12-18 symptoms. If the 
  number of obserations is large in both hospital and community samples,
  for exmaple, over 1000 cases total, use more symptoms, 
  otherwise use fewer. Sentivitity analysis can also 
  be used to choose \code{nsymp}. In general, the results stabilize in 
  the right range of the choices of \code{nsymp}. default=16.}
  \item{n.subset}{A positive integer specifing the total number of
     subsets and thus estimations of all symptoms.
     The default is \code{300}.}
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
   \item{boot.se}{a Logical value. If \code{TRUE}, bootstrap
     standard errors of the CSMF are estimated.  This typically takes a lot
    of computing time.  Default=\code{FALSE}.}
  \item{nboot}{a positive integer. If \code{boot.se=TRUE}, it
    specifies the number of bootstrapping samples taken to estimate the
    standard errors of CSMF. The default is \code{300}.}
  \item{printit}{Logical value. If \code{TRUE}, the progress of the
    estimation procedure will be printed on the screen.}
  \item{print.reg.size}{Logical value. If \code{TRUE}, the size of the
  regression matrix is printed at each step of subsampling. It provides
  helpful information for user to choose the number of symptoms to
  subsample. It is recommended to print the size of the regression
  matrix for different values of \texttt{nsymp} with a small size of
  \texttt{n.subset}. }
\item{predict.S}{Logical value. If \code{TRUE}, the predicted
  probabilities of each symptom in the community sample are estimated.
  If \code{boot.se} is \code{TRUE}, \code{predict.S} is a matrix with
  \code{nboot} rows and as many columns as the number of total symptoms
  used in \code{formula}.If \code{boot.se} is \code{FALSE},
  \code{predict.S} is a vector of the number of the symptoms. The
  default is \code{TRUE}.}
}
\details{
  For details, please refer to "Verbal Autposy Methods with
  Multiple Causes of Death"(King and Lu, 2008), and
  \url{http:\\gking.harvard.edu\va}
}

\value{
  \code{va} outputs a list containing the estimated cause-specific mortality
  fractions \code{est.CSMF}, and the true cause-specific mortality
  fractions \code{true.CSMF}, whenever available.

  If \code{boot.se=TRUE}, the bootstrapping estimations of
  \code{est.CSMF} and their standard errors \code{CSMF.se} are reported.

  When the causes of death are observed in validation studies,
  the bootstrapping mean(\code{true.CSMF.bootmean}) and
  standard error(\code{true.bootse}) of the sample CSMF are also reported.
}
\references{
  King, Gary and Ying Lu. (2008) \dQuote{Verbal Autopsy Methods with
    Multiple Causes of Death}, Statistical Science, 14(1).  Also
   available at \url{http:gking.harvard.edu/va}
}

\keyword{verbal autopsy, CSMF}
  
