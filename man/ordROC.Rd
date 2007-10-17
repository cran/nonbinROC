\name{ordROC}
\alias{ordROC}
\title{Performs ROC-type analysis for ordinal scale gold standard}
\description{
Assess the accuracies of diagnostic tests when the gold standard is ordinal.
}
\usage{
ordROC(gldstd, test1, test2 = NULL, penalty = NULL)
}
\arguments{
  \item{gldstd}{vector giving the gold standard}
  \item{test1}{vector giving the diagnostic test}
  \item{test2}{optional vector giving the diagnostic test}
  \item{penalty}{optional square matrix for the penalty function L[i,j] in which 0 <= L[i,j] <= 1 for
  j > i and the remaining elements can be set at any value. By default, L[i,j] = 1 for j > i and 0 otherwise.}
}

\value{
In the paired design, `ordROC' returns the pairwise accuracies and associated standard errors between 
each category of the gold standard, the penalty matrix, the overall accuracies and associated standard 
errors of the diagnostic tests, their overall covariance and the test for comparing two diagnostic tests' 
accuracies. 

Otherwise, `ordROC' returns the pairwise accuracies and associated standard errors between 
each category of the gold standard, the penalty matrix, the overall accuracy and associated standard 
error of the diagnostic test.
}

\references{
Obuchowski, N. A. (2005) Estimating and comparing diagnostic tests' accuracy when the gold 
standard is not binary. \emph{Academic Radiology}, 12, 1198-1204.}
\author{Paul Nguyen}

\examples{
#An example for heart tissue scarring from Obuchowski (2005)
data(heart)
attach(heart)
penalty = matrix(c(0,0,0,0,0.25,0,0,0,0.5,0.25,0,0,1,0.5,0.25,0), nrow = 4)
ordROC(PET, MRI, penalty = penalty)
}
\keyword{ univar }% at least one, from doc/KEYWORDS
