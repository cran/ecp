
// includes from the plugin

#include <Rcpp.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes


// declarations
extern "C" {
SEXP getWithin( SEXP alpha_, SEXP X_) ;
SEXP getBetween( SEXP alpha_, SEXP X_, SEXP Y_) ;
SEXP splitPointC( SEXP s_, SEXP e_, SEXP D_, SEXP min_size_) ;
}

// definition

SEXP getWithin( SEXP alpha_, SEXP X_ ){
BEGIN_RCPP
try{
NumericMatrix X(X_);//matrix of points used
		double alpha=as<double>(alpha_);
		double ret = 0.0;
		int n = X.nrow();
		for(int i=0;i<n;++i)
			for(int j=0;j<n;++j)
				ret += std::pow((sqrt(sum((X(i,_)-X(j,_))*(X(i,_)-X(j,_))))), alpha);
		return wrap(ret/(n*n));		
	}
catch(std::exception& ex){
		forward_exception_to_r(ex);
	}
	catch(...){
		Rf_error("unknown C++ exception");
	}
END_RCPP
}


SEXP getBetween( SEXP alpha_, SEXP X_, SEXP Y_ ){
BEGIN_RCPP

try{
		NumericMatrix X(X_), Y(Y_);
		double alpha=as<double>(alpha_);
		double ret = 0.0;
		int n = X.nrow(), m = Y.nrow();
		for(int i=0;i<n;++i)
			for(int j=0;j<m;++j)
				ret += std::pow((sqrt(sum((X(i,_)-Y(j,_))*(X(i,_)-Y(j,_))))), alpha);				
		return wrap(2*ret/(n*m));	
	}
	catch(std::exception& ex){
		forward_exception_to_r(ex);
	}
	catch(...){
		Rf_error("unknown C++ exception");
	}

END_RCPP
}


SEXP splitPointC( SEXP s_, SEXP e_, SEXP D_, SEXP min_size_ ){
BEGIN_RCPP

//This impementation takes at most O(n^2) time to find each change point.
//Thus if k change points are found, our algorithm is O(kn^2).
using namespace Rcpp;
//used to sum an individual row/colum of a matrix
#define SUM(A) std::accumulate(A.begin() , A.end() , 0.0)//uses STL templated funciton so need to use cxxfunciton()

NumericVector best = NumericVector::create(-1.0,R_NegInf);
int e = as<int>(e_), s = as<int>(s_), min_size = as<int>(min_size_);//convert to C++ primatives

NumericMatrix D(D_);//the distance matrix
e = e-s+1;//now represents the number of data points

int t1=min_size, t2=min_size<<1;//tau1 and tau2
NumericMatrix cut1 = D(Range(0,t1-1),Range(0,t1-1));
NumericMatrix cut2 = D(Range(t1,t2-1),Range(t1,t2-1));
NumericMatrix cut3 = D(Range(0,t1-1),Range(t1,t2-1));

double A = SUM(cut1)/2, B1=SUM(cut2)/2, AB1=SUM(cut3);
double tmp= 2*AB1/((t2-t1)*(t1)) - 2*B1/((t2-t1-1)*(t2-t1)) - 2*A/((t1-1)*(t1));
tmp *= (t1*(t2-t1)/t2);
if(tmp > best[1]){//update if better location is found
	best[0] = t1+s;
	best[1] = tmp;
}

t2+=1;

NumericVector B(e+1,B1), AB(e+1,AB1);//might be wasting a little space because of
//elements at the beginning that are not used
for(;t2<=e;++t2){//update between and within (for the right sample) distances
	B[t2] = B[t2-1] + SUM(D(Range(t2-1,t2-1),Range(t1,t2-2)));
	AB[t2] = AB[t2-1] + SUM(D(Range(t2-1,t2-1),Range(0,t1-1)));
	tmp = 2*AB[t2]/((t2-t1)*(t1))-2*B[t2]/((t2-t1-1)*(t2-t1))-2*A/((t1)*(t1-1));
	tmp *= (t1*(t2-t1)/t2);
	if(tmp > best[1]){//update if better location is found
		best[0] = t1+s;
		best[1] = tmp;
	}
}

t1+=1;

for(;;++t1){//iterate over remaining possible locations for tau1,tau2
	t2=t1+min_size;
	if(t2>e)//remaining interval is too small to fit another cluster
		break;
	double addA = SUM(D(Range(t1-1,t1-1),Range(0,t1-2)));
	A+=addA;//update within distance for left sample
	double addB = SUM(D(Range(t1-1,t1-1),Range(t1,t2-2)));
	for(;t2<=e;++t2){
		addB += D(t1-1,t2-1);
		B[t2]-=addB;
		AB[t2]+=(addB-addA);
		tmp = 2*AB[t2]/((t2-t1)*(t1))-2*B[t2]/((t2-t1-1)*(t2-t1))-2*A/((t1-1)*(t1));
		tmp *= (t1*(t2-t1)/t2);
		if(tmp > best[1]){
			best[0] = t1+s;
			best[1] = tmp;
		}
	}
}
return wrap(best);

END_RCPP
}



