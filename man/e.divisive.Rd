\name{e.divisive}
\alias{e.divisive}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	ENERGY DIVISIVE
}
\description{
	Divisive hierarchical algorithm for change point analysis.
}
\usage{
e.divisive(X, sig_lvl=.05, R=199, k=NULL, min_size=30, alpha=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  	A T by d matrix contating the time series with observations in R^d
	}
	\item{sig_lvl}{
		The level at which to test if a proposed change point is 
		statistically significant.
	}
	\item{R}{
		Number of permutations to use in the permutation test.
	}
	\item{k}{
		Number of change point locations to estimate. If this is NULL then only 
		the statistically significant estimated change points are returned.
	}
	\item{min_size}{
		Minimum number of observations between change points.
	}
	\item{alpha}{
		The index used for determining the distance between and within 
		segments.
	}
}
\details{
Segments are found through the use of a binary bisection method and a permutation 
test. The time complexity of this method is \emph{O(k*T^2)}, where \emph{k} is the 
number of change points found, and \emph{T} is the number of data points.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	The returned value is a list with the following components.
	\item{cluster_number}{The number of segments within the data created 
		by the change points.}
	\item{order_found}{The order in which the change points were found.}
	\item{list}{Locatoins of the statistically significant change points.}
	\item{considered_last}{Location of the last tested change point that 
		was not found to be statistically significant.}
	\item{pvals}{Approximate p-values found by each permutation test.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
James NA, Matteson DS (2013). A Nonparametric Approach for Multiple Change Point Analysis of Multivariate Data.

James NA, Matteson DS (2013). ecp: An R Package for Nonparametric Multiple Change Point Analysis of Multivariate Data.

Rizzo ML, Szekely GL (2005). Hierarchical clustering via joint between-within distances: Extending ward's minimum variance method. Journal of Classification. pp. 151 - 183.

Rizzo ML, Szekely GL (2010). Disco analysis: A nonparametric extension of analysis of variance. The Annals of Applied Statistics. pp. 1034 - 1055.
}
\author{
Nicholas A. James
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{e.agglo}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
set.seed(100)
x = matrix(c(rnorm(100),rnorm(100,3),rnorm(100,0,2)))
y = e.divisive(x)
y$list
# 1 101 201 301
y$cluster_number
# 3
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ divisive }
\keyword{ hierarchical }% __ONLY ONE__ keyword per line
