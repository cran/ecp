\name{e.cp3o}
\alias{e.cp3o}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	CHANGE POINTS ESTIMATION BY PROBABILISTICALLY PRUNED OBJECTIVE
}
\description{
	An algorithm for multiple change point analysis that uses dynamic programming and 
	probabilistic pruning.
}
\usage{
e.cp3o(Z, K=1, delta=29, alpha=1, eps=0.01, verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{Z}{
  		A T x d matrix containing the length T time series with d-dimensional observations.
  	}
	\item{K}{
		The maximum number of change points.
	}
	\item{delta}{
		The window size used to calculate the calculate the complete portion off our 
		approximate test statistic. This also corresponds to one less than the minimum 
		segment size.
	}
	\item{alpha}{
		The moment index used for determining the distance between and within 
		segments.
	}
	\item{eps}{
		The epsilon probability used for the probabilistic pruning procedure.
	}
	\item{verbose}{
		A flag indicating if status updates should be printed.
	}
}
\details{
Segmentations are found through the use of dynamic programming and probabilistic 
pruning. The computational complexity of this method is \emph{O(KT^2)}, where 
\emph{K} is the maximum number of change points, and \emph{T} is the number of 
observations.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	The returned value is a list with the following components.
	\item{number}{The estimated number of change points.}
	\item{estimates}{The location of the change points estimated by the procedure.}
	\item{gofM}{A vector of goodness of fit values for differing number of 
		change points. The first entry corresponds to when there is only 
		a single change point, the second for when there are two, and 
		so on.}
	\item{cpLoc}{A list of all the optimal change point locations for differing 
		numbers of change points. The first component corresponds to when there 
		is only one change point, the second for when there are two change points, 
		and so on.}
	\item{time}{The total amount to time take to estimate the change point locations.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~


Rizzo M.L., Szekely G.L (2005). Hierarchical clustering via joint between-within distances: Extending ward's minimum variance method. Journal of Classification.

Rizzo M.L., Szekely G.L. (2010). Disco analysis: A nonparametric extension of analysis of variance. The Annals of Applied Statistics.
}
\author{
Nicholas A. James
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
set.seed(400)
x1 = matrix(c(rnorm(100),rnorm(100,3),rnorm(100,0,2)))
y1 = e.cp3o(Z=x1, K=7, delta=29, alpha=1, eps=0.01, verbose=FALSE)
#View estimated change point locations
y1$estimates
x2 = rbind(MASS::mvrnorm(100,c(0,0),diag(2)),MASS::mvrnorm(100,c(2,2),diag(2)))
y2 = e.cp3o(Z=x2, K=4, delta=29, alpha=1, eps=0.01, verbose=FALSE)
#View estimated change point locations
y2$estimates
#View all possible segmentations for differing numbers of change points
y2$cpLoc
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pruning }
\keyword{ probabilistic }
\keyword{ dynamic }% __ONLY ONE__ keyword per line
