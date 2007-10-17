\name{nomROC}
\alias{nomROC}
\title{Performs ROC-type analysis for nominal scale gold standard}
\description{
Assess the accuracies of diagnostic tests when the gold standard is nominal.
}
\usage{
nomROC(gldstd, test1, test2 = NULL, penalty = NULL)
}

\arguments{
  \item{gldstd}{vector giving the gold standard}
  \item{test1}{data frame giving the set of confidence scores}
  \item{test2}{optional data frame giving another set of confidence scores}
  \item{penalty}{optional square matrix for the penalty function L[i,j] in which 0 <= L[i,j] <= 1 for
  j > i and the remaining elements can be set at any value. By default, L[i,j] = 1 for j > i and 0 otherwise.}
}
\value{
In the paired design, `nomROC' returns the pairwise accuracies and associated standard errors between 
each category of the gold standard, the penalty matrix, the overall accuracies and associated standard 
errors of the diagnostic tests, their overall covariance and the test for comparing two diagnostic tests' 
accuracies. 

Otherwise, `nomROC' returns the pairwise accuracies and associated standard errors between 
each category of the gold standard, the penalty matrix, the overall accuracy and associated standard 
error of the diagnostic test.
}

\references{
Obuchowski, N. A. (2005) Estimating and comparing diagnostic tests' accuracy when the gold 
standard is not binary. \emph{Academic Radiology}, 12, 1198-1204.}
\author{Paul Nguyen}

\examples{
#A paired design example for acute abdominal pain from Obuchowski, Goske and Applegate (2001)
data(abpain)
attach(abpain)
pre = data.frame(Pre1, Pre2, Pre3)
post = data.frame(Post1, Post2, Post3)
penalty = matrix(c(0,0,0,1,0,0,1,0.5,0), nrow = 3)
nomROC(Group, pre, post, penalty)
}

\keyword{ univar }% at least one, from doc/KEYWORDS
