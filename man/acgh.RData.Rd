\name{acgh}
\alias{acgh}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Bladder Tumor Micro-Array Data
}
\description{
Cleaned micro-array data for 43 different indivaiduals with a bladdar tumor.

}
\usage{data(acgh)}
\format{
  A list with the following components.
  
  data: The cleaned micro-array data for the 43 individuals. This information is stored in a 2215 by 43 matrix.
  
  individual: A numeric vector indicating which individuals' mico-array data remained after the cleaning process.
  
}
\source{
Bleakley K., Vert J.-P. (2011), The group fused Lasso for multiple change-point detection
}
\references{
Bleakley K., Vert J.-P. (2011), The group fused Lasso for multiple change-point detection

James NA, Matteson DS (2013). ecp: An R Package for Nonparametric Multiple Change Point Analysis of Multivariate Data.
}
\examples{
data(ecp/data/acgh.RData)
}
\keyword{datasets}
