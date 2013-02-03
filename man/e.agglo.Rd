\name{e.agglo}
\alias{e.agglo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	ENERGY AGGLOMERATIVE
}
\description{
	Agglomerative hierarchical algorithm for change point analysis.
}
\usage{
e.agglo(member, data, alpha=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{member}{Segment membership vector for the time series.}
  \item{data}{A matrix containing the time series with observations 
  	in R^d.}
  \item{alpha}{Index used for determining the distance between and
  	within segments.}
}
\details{
Segments are created based on the initial segmentation provided by the \emph{member} 
argument. In each iteration, segments are merged so as to maximize a goodness of fit 
statistic. The time complexity of this method is \emph{O(T^2)}, where \emph{T} is the 
number of data points.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list with the following components.
\item{merged}{A matrix indicating which segments were merged at each 
	step of the agglomerative procedure.}
\item{gof}{Vector showing the progression of the goodness of fit 
	statistic.}
\item{list}{Matrix showing the progression of the set of change points.}
\item{opt}{The location of the estimated change points.}
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
\code{\link{e.divisive}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
set.seed(100)
mem = rep(c(1,2,3,4),c(10,10,10,10))
x = c(rnorm(10),rnorm(20,2),rnorm(10,-1))
y = e.agglo(mem,as.matrix(x),1)
y$opt
#11 31
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ agglomerative }
\keyword{ hierarchical }% __ONLY ONE__ keyword per line
