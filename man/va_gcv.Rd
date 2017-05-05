\name{va.gcv}
\alias{va.gcv}

\title{va.gcv:General Cross Validation for Seleting Optimal Size to Subset Symptoms}
\description{
  This function determines the optimal size of the subsets of symptoms.
In general, the training set is randomy splited, half of the training
set is treated as the hospital sample and the other half is treated as
the community sample. \code{va.gcv} seraches for the optimal size of subsets
\code{nsymp} that minimize the prediction errors in the community
sample.The estimation is done using constrained quadratic optimization.}

\usage{

va.gcv(formula, data, nsymp.vec, n.subset=300, prob.wt=1,
  boot.se=FALSE, nboot=1, printit=FALSE, print.reg.size=TRUE)
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
  \item{nsymp.vec}{a vector of positive integer, it contains a collection 
    of different \code{nsymp} that can be used by \code{va()}. For
    a total of \code{J} number of causes of death and a total of
    \code{ns} symptoms in the sample, \code{nsymp.vec} cna be set to be
    a vector \code{a:b}, while \code{a} is the smallest integer than
    $2^a>J$. \code{b} is typically set to be \code{floor{0.75*b}}. If
    sample size is small, \code{b} can be set to smaller value to avoid
    function exiting due to data sparsity. No default value is set.} 
  \item{n.subset}{A positive integer specifing the total number of
     subsets and thus estimations of all symptoms.
     The default is \code{300}.}
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
     of computing time. It is highly suggested to set
     \code{boot.se=FALSE} in \code{va_gcv}. Default=\code{FALSE}.}
  \item{nboot}{a positive integer. If \code{boot.se=TRUE}, it
    specifies the number of bootstrapping samples taken to estimate the
    standard errors of CSMF. The default is \code{1}.}
  \item{printit}{Logical value. If \code{TRUE}, the progress of the
    estimation procedure will be printed on the screen.}
  \item{print.reg.size}{Logical value. If \code{TRUE}, the size of the
  regression matrix is printed at each step of subsampling. It provides
  helpful information for user to choose the number of symptoms to
  subsample. It is recommended to print the size of the regression
  matrix for different values of \texttt{nsymp} with a small size of
  \texttt{n.subset}. }
}
\details{
  For details, please refer to "Verbal Autposy Methods with
  Multiple Causes of Death"(King and Lu, 2008), and
  \url{http:\\gking.harvard.edu\va}
}

\value{
  \code{va.gcv} outputs two objects. \code{best.symp} returns the best
  \code{nsymp} that minimizes mean square error between estimated
  cause-specific mortality fraction and the observed cause-specific
  mortality fraction. \code{mse} returns a vector of mean square errors
  associated with each size of the subsets (as specified in \code{symp.vec}).
}
\references{
  King, Gary and Ying Lu. (2008) \dQuote{Verbal Autopsy Methods with
    Multiple Causes of Death}, 14(1), Statistical Science,  Also
   available at \url{http:gking.harvard.edu/va}
}

\keyword{verbal autopsy, CSMF}
  
