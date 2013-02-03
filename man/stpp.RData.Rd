\name{stpp}
\alias{stpp}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Spatio-temporal Point Process
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
A spatio-temporal inhomogeneous poisson point process. The data set contains 3 change points 
which occur at times 1, 3, and 4.5.
}
\usage{data(stpp)}
\format{
  A list with the following components.
  
  period1: The inhomogeneous poisson point process over the time interval [0,1).
  
  period2: The inhomogeneous poisson point process over the time interval [1,3).
  
  period3: The inhomogeneous poisson point process over the time interval [3,4.5).
  
  period4: The inhomogeneous poisson point process over the time interval [4.5,7].
  
  pointProcess: The inhomogeneous poisson point process over the time interval [0,7].
  
  member: The segment membership array for when the process is divided into 50 segments each spanning an equal amount of time.
}

\references{
James NA, Matteson DS (2013). ecp: An R Package for Nonparametric Multiple Change Point Analysis of Multivariate Data.
}
\examples{
data(ecp/data/stpp.RData)
}
\keyword{datasets}
