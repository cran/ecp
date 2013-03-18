\name{STPP}
\alias{STPP}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Spatio-temporal Point Process
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
A spatio-temporal inhomogeneous Poisson point process. The data set contains 3 change points 
which occur at times 1, 3, and 4.5.
}
\usage{data(STPP)}
\format{
  A list with the following components.
  
  period1: the inhomogeneous Poisson point process over the time interval [0,1).
  
  period2: the inhomogeneous Poisson point process over the time interval [1,3).
  
  period3: the inhomogeneous Poisson point process over the time interval [3,4.5).
  
  period4: the inhomogeneous Poisson point process over the time interval [4.5,7].
  
  pointProcess: the inhomogeneous Poisson point process over the time interval [0,7].
  
  member: the segment membership array for when the process is divided into 50 segments each spanning an equal amount of time.
}

\references{
James N.A., Matteson D.S. (2013). ecp: An R Package for Nonparametric Multiple Change Point Analysis of Multivariate Data.
}
\examples{
data(ecp/data/STPP.RData)
}
\keyword{datasets}
