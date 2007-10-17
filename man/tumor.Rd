\name{tumor}
\alias{tumor}
\docType{data}
\title{Renal tumor mass size}
\description{
A data frame giving a series of measurements of the renal tumor mass size for 74 patients 
based on surgery, computed tomography and a fictitious test.
}
\usage{data(tumor)}
\format{
  This data frame contains the following 3 variables measured in cm:
  \describe{
    \item{\code{SURG}}{renal mass size measured by surgery, considered as the gold standard}
    \item{\code{CT}}{renal mass size measured by CT images}
    \item{\code{Fi}}{renal mass size measured by a fictitious test}
  }
}
\source{
Obuchowski, N. A. (2005) Estimating and comparing diagnostic tests' accuracy when the gold 
standard is not binary. \emph{Academic Radiology}, 12, 1198-1204.
}
\keyword{datasets}
