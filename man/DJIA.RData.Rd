\name{DJIA}
\alias{DJIA}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Dow Jones Industrial Average Index
}
\description{
The weekly log returns for the Dow Jones Industrial Average index from April 1990 to January 2012.
}
\usage{data(DJIA)}
\format{
  A list with the following components.
  
  dates: A character vector of dates associated with each observation in the returns series.
  
  returns: Weekly log returns from April 1990 to January 2012.
}
\source{
http://research.stlouisfed.org/fred2/series/DJIA/downloaddata
}
\references{
James NA, Matteson DS (2013). ecp: An R Package for Nonparametric Multiple Change Point Analysis of Multivariate Data.
}
\examples{
data(ecp/data/DJIA.RData)
}
\keyword{datasets}
