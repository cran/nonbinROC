\name{heart}
\alias{heart}
\docType{data}
\title{Heart damage diagnosis after myocardial infarction}
\description{
A data frame giving a series of ordinal measurements of the positron emission tomography
and magnetic resonance imaging scans for the heart tissues of 241 fictitious 
patients after myocardial infarctions.
}
\usage{data(heart)}
\format{
  This data frame contains the following 2 variables:
  \describe{
    \item{\code{PET}}{gold standard for evaluating the most damage part of the heart 
    (1 = normal, 2 = ischemic, 3 = hibernating, and 4 = necrotic)}
    \item{\code{MRI}}{amount of scarring in the most damaged part of the heart (0 = normal, 
    1 = 1-24\%, 2 = 25-49\%, 3 = 50-74\%, 4 = 75-99\%, and 5 = 100\% scarring)}
  }
}
\source{
Obuchowski, N. A. (2005) Estimating and comparing diagnostic tests' accuracy when the gold 
standard is not binary. \emph{Academic Radiology}, 12, 1198-1204.
}
\keyword{datasets}
