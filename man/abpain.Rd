\name{abpain}
\alias{abpain}
\docType{data}
\title{Paediatric abdominal pain study}
\description{
A data frame giving a series of confidence scores for the pre-imaging and 
post-imaging diagnoses of 60 patients suffering from acute abdominal pain. 
}
\usage{data(abpain)}
\format{  
    This data frame contains the following 7 variables:
\describe{
    \item{\code{Group}}{gold standard for diagnosing the patient with acute abdominal pain
    (1 = surgical abdominal or urogenital condition, 2 = non-surgical abdominal condition and 
    3 = non-surgical urogenital condition)}
    \item{\code{Pre1}}{pre-image confidence score for surgical abdominal or urogenital condition}
    \item{\code{Pre2}}{pre-image confidence score for non-surgical abdominal condition}
    \item{\code{Pre3}}{pre-image confidence score for non-surgical urogenital condition}
    \item{\code{Post1}}{post-image confidence score for surgical abdominal or urogenital condition}
    \item{\code{Post2}}{post-image confidence score for non-surgical abdominal condition}
    \item{\code{Post3}}{post-image confidence score for non-surgical urogenital condition}
  }
}

\source{
Obuchowski, N. A., Goske M. J. and Applegate, K. E. (2001) Assessing physicians' accuracy in 
diagnosing paediatric patients with acute abdominal pain: measuring accuracy
for multiple diseases. \emph{Statistics in Medicine}, 20, 3261-3278.
}
\references{
Obuchowski, N. A. (2005) Estimating and comparing diagnostic tests' accuracy when the gold 
standard is not binary. \emph{Academic Radiology}, 12, 1198-1204.
}
\keyword{datasets}
