\name{contROC}
\alias{contROC}
\title{Performs ROC-type analysis for continuous scale gold standard}
\description{
Assess the accuracies of diagnostic tests when the gold standard is continuous.
}
\usage{
contROC(gldstd, test1, test2 = NULL)
}
\arguments{
  \item{gldstd}{vector giving the gold standard}
  \item{test1}{vector giving a diagnostic test}
  \item{test2}{optional vector giving another diagnostic test}
}
\value{
In the paired design, `contROC' returns the accuracies and associated standard errors of the 
diagnostic tests, their covariance and the test for comparing two diagnostic tests' accuracies. 

Otherwise, `contROC' returns the accuracy and associated standard error of the diagnostic test.
}
\references{
Obuchowski, N. A. (2005) Estimating and comparing diagnostic tests' accuracy when the gold 
standard is not binary. \emph{Academic Radiology}, 12, 1198-1204.}
\author{Paul Nguyen}

\examples{
#A paired design example for tumor renal masses from Obuchowski (2005)
data(tumor)
attach(tumor)
contROC(SURG, CT, Fi)

#A paired design example for blood iron concentrations from Obuchowski (2006)
data(blood)
attach(blood)
contROC(Serum, TIBC, TS)
}
\keyword{ univar }% at least one, from doc/KEYWORDS
