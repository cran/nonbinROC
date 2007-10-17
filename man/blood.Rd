\name{blood}
\alias{blood}
\docType{data}
\title{Blood iron concentration}
\description{
A data frame giving a series of measurements for the blood iron concentration of 55 anaemia female patients 
based on the serum ferritin concentration, percent transferrin saturation and total iron binding capacity. 
}
\usage{data(blood)}
\format{
  This data frame contains the following 3 variables:
  \describe{
    \item{\code{Serum}}{blood iron concentration measured by the serum ferritin concentration, 
        considered as the gold standard}
    \item{\code{TIBC}}{blood iron concentration measured by total iron binding capacity}
    \item{\code{TS}}{blood iron concentration measured by percent transferrin saturation}
  }
}

\source{
Obuchowski, N. A. (2006) An ROC-type measure of diagnostic accuracy when the gold standard is continuous-scale. 
\emph{Statistics in Medicine}, 25, 481-493.
}
\keyword{datasets}
