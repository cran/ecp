#include<Rcpp.h>
#include<math.h>

#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

// user includes
#include<algorithm>
#include<iostream>
#include<vector>
#include<set>
// #include<random>

using namespace std;
using namespace Rcpp;

// declarations
extern "C"{
SEXP eFastC(SEXP Z_, SEXP K_, SEXP delta_, SEXP alpha_, SEXP eps_, SEXP verbose_);
SEXP ksFastC(SEXP Z_, SEXP K_, SEXP minsize_, SEXP eps_, SEXP verbose_);
}

double dst(const NumericVector& X, const NumericVector& Y, double alpha);
std::vector<double> MEAN_VAR(const std::vector<double>& X);
int MSample(int a, int b);
double get_gamma(int N, double eps, double minsize, const NumericVector& DLL,
                 const NumericVector& DRR, const NumericVector& DLR,
	         const NumericVector& cSum);
double delta_sum(NumericMatrix& X, int a, int b, double alpha);
std::vector<std::vector<int> > find_locations(const NumericMatrix& A);
double get_gamma_ks(int N, double eps, double minsize, NumericMatrix& Z);
double dist_ks(std::vector<double>& X, std::vector<double>& Y);

// definitions

SEXP eFastC( SEXP Z_, SEXP K_, SEXP delta_, SEXP alpha_, SEXP eps_, SEXP verbose_){
BEGIN_RCPP
	
	//Convert SEXP variables to types Rcpp/C++ knows
	int K = as<int>(K_), delta = as<int>(delta_);
	double alpha = as<double>(alpha_), eps = as<double>(eps_);
	bool verbose = as<bool>(verbose_);
	NumericMatrix Z(Z_);
	
	int N = Z.nrow(); //Get number of observations

	//Since the update of the goodness-of-fit value with k change points only depends
	//upon the goodness-of-fit values with (k-1) change points we can use a 2xN matrix.
	//However, we will need a KxN matrix to store the change point locaitons.
	NumericMatrix FF(2,N), A(K,N); //Create matrices used to store goodness-of-fit values and 
	//change point locations

	//Goodness-of-fit values for when outer loop is at the end of the time series.
	//GOF[k] corresponds to the goodness-of-fit value  for the entire time series when
	//using k+1 change points.
	NumericVector GOF(K); 

	std::fill(FF.begin(), FF.end(), R_NegInf); //Fill FF matrix with negative infinity
	std::fill(A.begin(), A.end(), -1); //Fill A matrix with -1 because of C arrays being 0 indexed

	//Counters and distances

	//Sum of within sample for left and right segments as well as between segment distances
	double dll = 0.0, drr = 0.0, dlr = 0.0;
	//Counter for number of terms in each sum
	int cll = 0, crr = 0, clr = 0;

	//Precalculations
	if(verbose)
		Rcout<<"Starting pre-processing"<<std::endl;
	int minsize = delta + 1;
	
	//Calculations for complete portion of statistics in delta neighborhoods
	NumericVector DLL(N), DRR(N), DLR(N);
	//DRR[s] = within sum for Z[s+1], Z[s+2], ..., Z[s+delta]
	//DLL[s] = within sum for Z[s-delta+1], Z[s-delta+2], ..., Z[s]
	//DLR[s] = between sum for the sets used to calculate DLL[s] and DRR[s]
	for(int s = delta; s < N; ++s){
		DLL[s] = delta_sum(Z, s-delta+1, s, alpha);
		if(s >= delta)//avoid array out of bounds error
			DRR[s-delta] = DLL[s];
	}

	//Calculate DLR array in O(delta*N) time

	NumericMatrix Left(N,2), Right(N,2);
	//Left(i,0) = sum of distances of Z[i] to {Z[i-1], Z[i-2], ..., Z[i-delta]}
	//Right(i,0) = sum of distances of Z[i] to {Z[i+1], Z[i+2], ..., Z[i+delta]}
	//Left(i,1) = sum of distances of Z[i] to {Z[i-delta], Z[i-delta-1], ..., Z[i-2*delta+1]}
	//Right(i,1) = sum of distances of Z[i] to {Z[i+delta], Z[i+delta+1], ..., Z[i+2*delta-1]}

	for(int i = delta; i < N-delta; ++i){
		for(int j1 = i-delta; j1 < i; ++j1)
			Left(i,0) += dst(Z(i,_), Z(j1,_), alpha);
		for(int j2 = i+1; j2 < i+delta+1; ++j2)
			Right(i,0) += dst(Z(i,_), Z(j2,_), alpha);

		//Avoid array out of bounds errors
		if(i >= 2*delta-1)
			for(int j3=i-2*delta+1; j3<i-delta+1; ++j3)
				Left(i,1) += dst(Z(i,_), Z(j3,_), alpha);
		if(i+2*delta-1 < N)
			for(int j4=i+delta; j4<i+2*delta; ++j4)
				Right(i,1) += dst(Z(i,_), Z(j4,_), alpha);
	}

	//Update DLR
	for(int i = 1; i < minsize; ++i)
		for(int j = minsize; j < minsize+delta; ++j)
			DLR[minsize-1] = DLR[minsize-1] + dst(Z(i,_), Z(j,_), alpha);
	
	for(int s = minsize; s < N-delta; ++s){
		double r1 = Left(s,0); //Z[s] moved from the right to left segment so track affected distances
		double r2 = Right(s-delta,1); //Z[s-delta] has been removed from the sample under consideration
		double r3 = dst(Z(s,_),Z(s+delta,_),alpha); //Account for double counting in a1 and a2 below
		double a1 = Left(s+delta,1); //Z[s+delta] has been added to the sample under consideration
		double a2 = Right(s,0); //Z[s] moved from the right to left segment so track affected distances
		double a3 = dst(Z(s,_),Z(s-delta,_),alpha); //Account for double counting in r1 and r2 above
		DLR[s] = DLR[s-1] - r1 - r2 - r3 + a1 + a2 + a3;
	}

	//Calculaton of cumulative sum of distances for adjacent observations
	NumericVector cSum(N);
	std::fill(cSum.begin(), cSum.end(), 0.0);
	//cSum[i] = |Z[0]-Z[1]|^alpha + |Z[1]-Z[2]|^alpha + ... + |Z[i-1]-Z[i]|^alpha
	for(int i = 1;i < N; ++i)
		cSum[i] = cSum[i-1] + dst(Z(i,_), Z(i-1,_), alpha);
	
	//Calculate the value of Gamma used for pruning
	double Gma = get_gamma(N, eps, minsize, DLL, DRR, DLR, cSum);

	if(verbose)
		Rcout<<"Pre-processing complete. Starting optimization."<<std::endl;

	//Solve for K=1
	std::vector<std::set<int> > testPoints(N);
	//testPoints[i] holds the pruned set of change point locations before i
	//Use std::set to make sure that there are no duplicate elements, plus have fast 
	//deletion of elements.
	for(int t = 2*minsize-1; t < N; ++t){//Iterate over "ending position" for time series
		std::set<int> cpSet; //Initially all locations are posible change point locations
		for(int a=delta; a<=t-minsize; ++a)
			cpSet.insert(a);
		for(std::set<int>::iterator si=cpSet.begin(); si!=cpSet.end(); ++si){
		//Iterate over possible change point locations for a given "ending position"
			int s = *si;//change point location under consideration
			int u = -1;//Only allowed to have 1 change point in this case (K=1)
			cll = crr = delta*(delta-1)/2;//Initialize counters
			clr = delta*delta;

			dll = DLL[s];//Initial within distnace sum for left segment
			drr = DRR[s];//Initial within distance sum for right segment
			dlr = DLR[s];//Initial between distance sum
			//Augment distnace sums
			dll += cSum[s-delta];
			drr += ( cSum[t] - cSum[s+delta] );
			cll += (s-delta-u);
			crr += (t-s-delta);

			//Calculate test statistic
			double stat = 2.0*dlr/(clr+0.0) - drr/(crr+0.0) - dll/(cll+0.0);
			double num = (s-u)*(t-s)*(1.0);
			double dnom = std::pow(t-u+0.0,2.0);
			stat *= (num/dnom);
			if(stat > FF(0,t)){//Update goodness-of-fit and location matrices
				FF(0,t) = stat;
				A(0,t) = s;
			}
		}

		//Update the set of possible change point locations for a given outer loop value
		std::set<int> cpSetCopy = cpSet;
		//Iterate over the copy and remove elements from the original
		for(std::set<int>::iterator si=cpSetCopy.begin(); si!=cpSetCopy.end(); ++si){
			int a = *si;
			int b = A(0,a);
			double V;
			if(b <= 0)
				V = 0.0;
			else
				V = cSum[b-1];
			double cst = (a-b)*(t-a)/std::pow(t-b+0.0,2.0);
			double stat = 2*DLR[a]/(delta*delta);
			double t1 = (DRR[a]+cSum[t]-cSum[a-delta])/(delta*(delta-1)/2+t-a-delta);
			double t2 = (DLL[a]+cSum[a-delta]-V)/(delta*(delta-1)/2+a-b-delta);
			stat -= (t1+t2);
			stat *= cst;
			//Check pruning condition and remove elements as necessary
			if(FF(0,a) + stat + Gma < FF(0,t))
				cpSet.erase(a);
		}
		cpSet.insert(delta);
		cpSet.insert(t-minsize);
		testPoints[t] = cpSet;
	}

	GOF[0] = FF(0,N-1);

	//If only 1 change point is to be found terminate early
	if(K == 1){
		if(verbose)
			Rcout<<"Finished optimization."<<std::endl;
		return wrap(List::create(_["number"]=1, _["estimates"]=A(0,N-1), _["gofM"]=GOF, _["cpLoc"]=A));
	}

	//Solve for K > 1 cases
	//Since FF only has 2 rows we will use the variable 'flag' to help switch between the 0th and 1st row
	bool flag = true;

	for(int k=1; k<K; ++k){
		//Find change points
		for(int t=2*minsize-1; t<N; ++t){
			FF(flag,t) = R_NegInf; //Reset goodness of fit value
			std::set<int> cpSet = testPoints[t];
			for(std::set<int>::iterator si=cpSet.begin(); si!=cpSet.end(); ++si){
				int s = *si;
				int u = A(k-1,s); //location of change point before s
				cll = crr = delta*(delta-1)/2;//initialize counters
				clr = delta*delta;

				dll = DLL[s]; //initial within distance sum for left segment
				drr = DRR[s]; //initial within distnace sum for right segment
				dlr = DLR[s]; //initial between distance sum

				//Augment distance sums
				dll += cSum[s-delta];
				if(u>0)
					dll -= cSum[u-1];
				drr += (cSum[t] - cSum[s+delta]);
				cll += (s-delta-u);
				crr += (t-s-delta);

				double stat = 2*dlr/(clr+0.0) - drr/(crr+0.0) - dll/(cll+0.0);
				double num = (s-u)*(t-s)*1.0;
				double dnom = std::pow(t-u+0.0, 2.0);
				stat *= (num/dnom);
				if(u>0) //Optimal partition cost if previous change points exists
					stat += FF(1-flag,s);
				if(stat > FF(flag,t)){
					FF(flag,t) = stat;
					A(k,t) = s;
				}
			}
			//Update the set of possible change point locations for a given set of outer loop values
			std::set<int> cpSetCopy = cpSet;
			for(std::set<int>::iterator si=cpSetCopy.begin(); si!=cpSetCopy.end(); ++si){
				int a = *si;
				int b = A(k,a);
				double V;
				if(b <= 0)
					V = 0;
				else
					V = cSum[b-1];
				double cst = (a-b)*(t-a)*1.0/((t-b)*(t-b));
				double stat = 2*DLR[a]/(delta*delta);
				double t1 = (DRR[a]+cSum[t]-cSum[a-delta])/(delta*(delta-1)/2+t-a-delta);
				double t2 = (DLL[a]+cSum[a-delta]-V)/(delta*(delta-1)/2+a-b-delta);
				stat -= (t1+t2);
				stat *= cst;
				if(FF(flag,a) + stat + Gma < FF(flag,t))
					cpSet.erase(a);
			}
			cpSet.insert(delta);
			cpSet.insert(t-minsize);
			testPoints[t] = cpSet;
		}

		GOF[k] = FF(flag,N-1); //Obtain best goodness-of-fit value for the entire series
		flag = 1-flag;//Flip value of flag so that it now points to the alternate row
	}

	if(verbose)
		Rcout<<"Finished Optimization."<<std::endl;

	//For some reason R/Rcpp doesn't like the use use std::adjacent_difference so 
	//had to explicitly do it myself.
	std::vector<double> f;
	for(int i=1; i<GOF.size(); ++i)
		f.push_back( GOF[i] - GOF[i-1] );
	//Calculate the sample mean and variance of vector f
	std::vector<double> mean_and_var = MEAN_VAR(f);
	double u = mean_and_var[0] + 0.5*std::pow( mean_and_var[1], 0.5 );
	int k = 0;
	//Check stopping rule to determine the number of change points
	for(int i=0; i<f.size() && f[i] >= u; ++i, ++k);
	
	//Find all optimal segmentations for differing numbers 
	//of change points.
	std::vector<std::vector<int> > cpLocs = find_locations(A);

	return wrap(List::create(_["number"]=k+1, _["estimates"]=cpLocs[k], _["gofM"]=GOF, _["cpLoc"]=cpLocs, _["gamma"]=Gma));

END_RCPP
}


std::vector<double> MEAN_VAR(const std::vector<double>& X){
	//Calculate the mean and sample variane for observations in X
	//The variance is calculated in a single pass
	if(X.size() == 0){
		std::vector<double> ret;
		ret.push_back(0); ret.push_back(0);
		return ret;
	}

	double mean = *(X.begin()), var = 0;

	std::vector<double>::const_iterator i = X.begin();
	++i;
	int cnt = 2;
	for(; i!=X.end(); ++i, ++cnt){
		double dif = *i-mean;
		var = var*(cnt-2)/((double)cnt-1) + dif*dif/cnt;
		mean += dif/cnt;
	}
	
	std::vector<double> res;
	res.push_back(mean); res.push_back(var);
	return res;
}

double dst(const NumericVector& X, const NumericVector& Y, double alpha){
	//Calculate Euclidean distance between X and Y
	NumericVector res = X-Y;
	double ip = std::inner_product(res.begin(), res.end(), res.begin(),0.0);
	return std::pow( ip, alpha/2 );
}


double get_gamma(int N, double eps, double minsize, const NumericVector& DLL,
                 const NumericVector& DRR, const NumericVector& DLR,
		 const NumericVector& cSum){
	//Calculate the value of gamma used in the pruning step
	double delta = minsize-1;
	//Determine the number of random samples to use
	int R = (int)(std::pow((double)N,0.5)/eps)+1;
	if(R < 100)//perform at least 100 samples
		R = 100;
	NumericVector kVec(R);
	int i=0;
	while(i<R){
		//Draw segment boundaries
		//v<t<s<u with t-v, s-t, and u-s all at least minsize 
		// Since the distribution for which such points can be uniformly
		// drawn from is complicated it is easier (and possibly faster)
		// to use acceptance/rejection.
		NumericVector samp(4);
		samp[0] = MSample(minsize, N-minsize) - 1;
		samp[1] = MSample(minsize, N-minsize) - 1;
		samp[2] = MSample(minsize, N-minsize) - 1;
		samp[3] = MSample(minsize, N-minsize) - 1;
		std::sort(samp.begin(), samp.end());
		int v = samp[0];
		int t = samp[1];
		int s = samp[2];
		int u = samp[3];
		//subtract 1 because C arrays are zero based
		int tv = t-v;
		int st = s-t;
		int us = u-s;
		if(tv>=minsize && st>=minsize && us>=minsize){
			double V;
			if(v <= 1)
				V = 0;
			else
				V = cSum[v-1];
			//Calculate the necessary differences
			double cst1 = (double)(t-v)*(u-t)/((u-v)*(u-v));
			double cst2 = (double)(t-v)*(s-t)/((s-v)*(s-v));
			double cst3 = (double)(s-t)*(u-s)/((u-t)*(u-t));

			double stat1 = 2*DLR[t]/(delta*delta) - (DRR[t]+cSum[u]-cSum[t+delta])/(delta*(delta-1)/2+u-t-delta) - (DLL[t]+cSum[t-delta]-V)/(delta*(delta-1)/2+t-v-delta);
			double stat2 = 2*DLR[t]/(delta*delta) - (DRR[t]+cSum[s]-cSum[t+delta])/(delta*(delta-1)/2+s-t-delta) - (DLL[t]+cSum[t-delta]-V)/(delta*(delta-1)/2+t-v-delta);
			double stat3 = 2*DLR[s]/(delta*delta) - (DRR[s]+cSum[u]-cSum[s+delta])/(delta*(delta-1)/2+u-s-delta) - (DLL[s]+cSum[s-delta]-cSum[t-1])/(delta*(delta-1)/2+s-t-delta);

			stat1 *= cst1; stat2 *= cst2; stat3 *= cst3;
			kVec[i] = stat1 - stat2 - stat3;
			i+=1;
		}
	}
	//Find 1-eps quantile estimate
	std::nth_element(kVec.begin(), kVec.begin() + (int)((1-eps)*R), kVec.end());
	return kVec[(int)((1-eps)*R)];
}

int MSample(int a, int b)
{
      GetRNGstate();
      // Generate a unifrom random number in (a, b+1).
      double r = a + (b + 1 - a) * unif_rand();
      PutRNGstate();
      // Take the integer part of the generated number.
      return (int)r;
}

double delta_sum(NumericMatrix& X, int a, int b, double alpha){
	//Determine the sum of the alpha distances between X[a], X[a+1], ..., X[b]
	double ret = 0.0;
	for(int i=a; i<b; ++i)
		for(int j=i+1; j<=b; ++j)
			ret += dst(X(i,_), X(j,_), alpha);
	
	return ret;
}

std::vector<std::vector<int> > find_locations(const NumericMatrix& A){
	//Determine all of the optimal segmentations for 
	//differing numbers of change points
	std::vector<std::vector<int> > res;
	//Obtain dimensions for matrix A
	//N = number of observations, K = maximum number of fit change points
	int K = A.nrow(), N = A.ncol();
	//k+1 is the number of change points
	for(int k=0; k<K; ++k){
		int cp = A(k,N-1), k1 = k;
		std::vector<int> cps;
		do{
			cps.push_back(cp+1);
			--k1;
			if(k1 >= 0)
				cp = A(k1,cp);
		} while(k1 >= 0);
		sort(cps.begin(),cps.end());
		res.push_back(cps);
	}
	return res;
}

SEXP ksFastC( SEXP Z_, SEXP K_, SEXP minsize_, SEXP eps_, SEXP verbose_){
BEGIN_RCPP
	
	//Convert SEXP variables to types Rcpp/C++ knows
	int K = as<int>(K_), minsize = as<int>(minsize_);
	double eps = as<double>(eps_);
	bool verbose = as<bool>(verbose_);
	NumericMatrix Z(Z_);
	
	int N = Z.nrow(); //Get number of observations

	//Since the update of the goodness-of-fit value with k change points only depends
	//upon the goodness-of-fit values with (k-1) change points we can use a 2xN matrix.
	//However, we will need a KxN matrix to store the change point locaitons.
	NumericMatrix FF(2,N), A(K,N); //Create matrices used to store goodness-of-fit values and 
	//change point locations

	//Goodness-of-fit values for when outer loop is at the end of the time series.
	//GOF[k] corresponds to the goodness-of-fit value  for the entire time series when
	//using k+1 change points.
	NumericVector GOF(K); 

	std::fill(FF.begin(), FF.end(), R_NegInf); //Fill FF matrix with negative infinity
	std::fill(A.begin(), A.end(), -1); //Fill A matrix with -1 because of C arrays being 0 indexed

	//Precalculations
	if(verbose)
		Rcout<<"Starting pre-processing"<<std::endl;
	
	//Calculate the value of Gamma used for pruning
	double Gma = get_gamma_ks(N, eps, minsize, Z);

	if(verbose)
		Rcout<<"Pre-processing complete. Starting optimization."<<std::endl;

	//Solve for K=1
	std::vector<std::set<int> > testPoints(N);
	//testPoints[i] holds the pruned set of change point locations before i
	//Use std::set to make sure that there are no duplicate elements, plus have fast 
	//deletion of elements.
	std::vector<std::vector<pair<int,double> > > ksVal;
	ksVal.reserve(N-2*minsize+1);
	// ksVals contains the test statistics to be used for the next iteration of K
	// ksVals[t] holds the test statistics corresponding to change points in textPoints[t]
	// initializing
	std::set<int> cpSet;
	std::vector<double> vec_us_sort(minsize);
	vec_us_sort.reserve(N-minsize);
	std::vector<pair<double,int> > vec_st_sort;
	vec_st_sort.reserve(N-minsize);
	for(int t=2*minsize-1; t<N; ++t){//Iterate over "ending position" for time series
		if(t==2*minsize-1){// initialize at first t
			//Initially all locations are posible change point locations
			for(int a=minsize-1; a<=t-minsize; ++a)
				cpSet.insert(a);
			// initializing samples at first change point
			// s=minsize-1, u=-1
			for(int i=0; i<=minsize-1; ++i)
				vec_us_sort.push_back(Z[i]);
			for(int i=minsize; i<=t; ++i)
				vec_st_sort.push_back(std::make_pair(Z[i],i));
			// sorting. vec_st need to be resorted at each t
			std::sort(vec_us_sort.begin(),vec_us_sort.end());
			std::sort(vec_st_sort.begin(),vec_st_sort.end());
		}else{ // each iteration, add one element to series under consideration
			cpSet.insert(t-minsize);
			vec_st_sort.push_back(std::make_pair(Z[t],t));
			std::sort(vec_st_sort.begin(),vec_st_sort.end());
		}
		// initializing vector before and after change point
		std::vector<double> vec_us = vec_us_sort; 
		std::vector<pair<double,int> > vec_st_pair = vec_st_sort;
		std::vector<double> vec_st;
		vec_st.reserve(N-minsize);
		// initializing vector containing test stat for pruning
		std::vector<double> stat_prune;
		stat_prune.reserve(cpSet.size());
		//Iterate over possible change point locations for a given "ending position" t
		for(std::set<int>::iterator si=cpSet.begin(); si!=cpSet.end(); ++si){
			int s = *si;//change point location under consideration
			int u = -1;//Only allowed to have 1 change point in this case (K=1)			
			if(s==minsize-1){// initialize at first change point
				// find the vector before and after change point
				for(int i=0; i<=vec_st_sort.size(); ++i){
					vec_st.push_back(vec_st_sort[i].first);
				}
			}else{// each iteration, change point location will proceed by one, so vector before and after change point need to be updated
				// find index to be removed
				for(int to_erase=0; to_erase<=vec_st_pair.size(); ++to_erase){
					if(vec_st_pair[to_erase].second==s){
						vec_us.push_back(vec_st_pair[to_erase].first);
						std::sort(vec_us.begin(),vec_us.end());
						vec_st_pair.erase(vec_st_pair.begin()+to_erase);
						vec_st.erase(vec_st.begin()+to_erase);
						break;
					}
				}
			}
			//Calculate test statistic		
			double stat = dist_ks(vec_us,vec_st);
			if(stat > FF(0,t)){//Update goodness-of-fit and location matrices
				FF(0,t) = stat;
				A(0,t) = s;
			}
			//Calculate test statistic for pruning
			if(s<2*minsize-1){
				stat_prune.push_back(stat);
			}else{
				int b = A(0,s);
				vector<double> vec_bs(s-b);
				for(int i=b+1; i<=s; ++i)
					vec_bs[i-b-1] = Z[i];
				std::sort(vec_bs.begin(),vec_bs.end());
				stat_prune.push_back(dist_ks(vec_bs,vec_st));
			}
		}
		//Update the set of possible change point locations for a given outer loop value
		std::set<int> cpSetCopy = cpSet;
		std::vector<pair<int,double> > ksVal_t; // element of ksVal
		ksVal_t.reserve(cpSet.size());
		//Iterate over the original and remove elements from the copy
		for(std::set<int>::iterator si=cpSet.begin(); si!=cpSet.end(); ++si){
			int a = *si;
			//Check pruning condition and remove elements as necessary
			if(FF(0,a) + stat_prune[a-minsize+1] + Gma < FF(0,t))
				cpSetCopy.erase(a);
			else
				ksVal_t.push_back(std::make_pair(a,stat_prune[a-minsize+1]));
		}
		ksVal.push_back(ksVal_t); // may not contain change points minsize-1 and t-minsize
		cpSetCopy.insert(minsize-1);
		cpSetCopy.insert(t-minsize);
		testPoints[t] = cpSetCopy;
	}

	GOF[0] = FF(0,N-1);

	//If only 1 change point is to be found terminate early
	if(K == 1){
		if(verbose)
			Rcout<<"Finished optimization."<<std::endl;
		return wrap(List::create(_["number"]=1, _["estimates"]=A(0,N-1), _["gofM"]=GOF, _["cpLoc"]=A));
	}

	//Solve for K > 1 cases
	//Since FF only has 2 rows we will use the variable 'flag' to help switch between the 0th and 1st row
	bool flag = true;

	//Initialize vector before change point minsize-1
	std::vector<double> vec_us1(minsize);
	for(int i=0; i<=minsize-1; ++i)
		vec_us1[i] = Z[i];
	std::sort(vec_us1.begin(),vec_us1.end());
	for(int k=1; k<K; ++k){
		std::vector<std::vector<pair<int,double> > > ksValNew;
		ksValNew.reserve(N-2*minsize+1);
		//Initialize vector after change point minsize-1
		std::vector<double> vec_st1;
		vec_st1.reserve(N-minsize);
		for(int i=minsize; i<=2*minsize-2; ++i)
			vec_st1.push_back(Z[i]);
		//Find change points
		for(int t=2*minsize-1; t<N; ++t){
			vec_st1.push_back(Z[t]);
			FF(flag,t) = R_NegInf; //Reset goodness of fit value
			std::set<int> cpSet = testPoints[t];
			// retrieve test stat calculated in (k-1)th iteration
			std::vector<pair<int,double> > ksVal_t = ksVal[t-2*minsize+1];
			// initializing vector containing test stat for pruning
			std::vector<double> stat_prune;
			stat_prune.reserve(cpSet.size());
			for(std::set<int>::iterator si=cpSet.begin(); si!=cpSet.end(); ++si){
				int s = *si;
				int u = A(k-1,s); //location of change point before s
				// retrieve test stat calculated in (k-1)th iteration
				int j=0; // find location of saved test stat
				double stat=0.0;				
				if(ksVal_t.size()>0){
					for(j=0; j<ksVal_t.size(); ++j){
						if(ksVal_t[j].first==s)
							break;
					}
				}		
				//Calculate test statistic for pruning
				if(ksVal_t.size()>0 && ksVal_t[j].first==s){
					stat = ksVal_t[j].second;
					ksVal_t.erase(ksVal_t.begin()+j);
					//Calculate test statistic for pruning
					int b = A(k,s);
					std::vector<double> vec_bs(s-b),vec_st(t-s);
					for(int i=b+1; i<=s; ++i)
						vec_bs[i-b-1] = Z[i];
					for(int i=s+1; i<=t; ++i)
						vec_st[i-s-1] = Z[i];
					std::sort(vec_bs.begin(),vec_bs.end());
					std::sort(vec_st.begin(),vec_st.end());
					stat_prune.push_back(dist_ks(vec_bs,vec_st));
				}else{
					if(s==minsize-1){
						std::sort(vec_st1.begin(),vec_st1.end());
						stat = dist_ks(vec_us1,vec_st1);
						//Calculate test statistic for pruning
						stat_prune.push_back(stat);
					}else{ // s==t-minsize
						std::vector<double> vec_us2(s-u), vec_st2(t-s);
						for(int i=u+1; i<=s; ++i)
							vec_us2[i-u-1] = Z[i];
						for(int i=s+1; i<=t; ++i)
							vec_st2[i-s-1] = Z[i];
						std::sort(vec_us2.begin(),vec_us2.end());
						std::sort(vec_st2.begin(),vec_st2.end());
						stat = dist_ks(vec_us2,vec_st2);
						//Calculate test statistic for pruning
						int b = A(k,s);
						std::vector<double> vec_bs(s-b);
						for(int i=b+1; i<=s; ++i)
							vec_bs[i-b-1] = Z[i];
						std::sort(vec_bs.begin(),vec_bs.end());
						stat_prune.push_back(dist_ks(vec_bs,vec_st2));
					}
				}
				if(u>0) //Optimal partition cost if previous change points exists
					stat += FF(1-flag,s);
				if(stat > FF(flag,t)){
					FF(flag,t) = stat;
					A(k,t) = s;
				}
			}
			//Update the set of possible change point locations for a given set of outer loop values
			std::set<int> cpSetCopy = cpSet;
			std::vector<pair<int,double> > ksValNew_t;
			ksValNew_t.reserve(cpSet.size());
			int prune_counter = 0;
			for(std::set<int>::iterator si=cpSetCopy.begin(); si!=cpSetCopy.end(); ++si){
				int a = *si;
				if(FF(flag,a) + stat_prune[prune_counter] + Gma < FF(flag,t))
					cpSet.erase(a);
				else
					ksValNew_t.push_back(std::make_pair(a,stat_prune[prune_counter]));	
				++prune_counter;			
			}
			ksValNew.push_back(ksValNew_t); // may not contain change points minsize-1 and t-minsize
			cpSet.insert(minsize-1);
			cpSet.insert(t-minsize);
			testPoints[t] = cpSet;
		}
		std::vector<std::vector<pair<int,double> > > ksVal = ksValNew;

		GOF[k] = FF(flag,N-1); //Obtain best goodness-of-fit value for the entire series
		flag = 1-flag;//Flip value of flag so that it now points to the alternate row
	}

	if(verbose)
		Rcout<<"Finished Optimization."<<std::endl;

	//For some reason R/Rcpp doesn't like the use use std::adjacent_difference so 
	//had to explicitly do it myself.
	std::vector<double> f;
	for(int i=1; i<GOF.size(); ++i)
		f.push_back( GOF[i] - GOF[i-1] );
	//Calculate the sample mean and variance of vector f
	std::vector<double> mean_and_var = MEAN_VAR(f);
	double u = mean_and_var[0] + 0.5*std::pow( mean_and_var[1], 0.5 );
	int k = 0;
	//Check stopping rule to determine the number of change points
	for(int i=0; i<f.size() && f[i] >= u; ++i, ++k);
	
	//Find all optimal segmentations for differing numbers 
	//of change points.
	std::vector<std::vector<int> > cpLocs = find_locations(A);

	return wrap(List::create(_["number"]=k+1, _["estimates"]=cpLocs[k], _["gofM"]=GOF, _["cpLoc"]=cpLocs, _["gamma"]=Gma));

END_RCPP
}


double get_gamma_ks(int N, double eps, double minsize, NumericMatrix& Z){
	//Determine the number of random samples to use
	int R = (int)(std::pow((double)N, 0.5)/eps)+1;
	if(R < 100)//perform at least 100 samples
		R = 100;
	NumericVector kVec(R);
	int j=0;
	while(j<R){
		// Draw segment boundaries
		// v<t<s<u with t-v, s-t, and u-s all at least minsize.
		// Since the distribution for which such points can be uniformly
		// drawn from is complicated it is easier (and possibly faster)
		// to use acceptance/rejection.
		NumericVector samp(4);
		//subtract 1 because C arrays are zero based
		samp[0] = MSample(minsize, N-minsize) - 1;
		samp[1] = MSample(minsize, N-minsize) - 1;
		samp[2] = MSample(minsize, N-minsize) - 1;
		samp[3] = MSample(minsize, N-minsize) - 1;
		std::sort(samp.begin(), samp.end());
		int v = samp[0];
		int t = samp[1];
		int s = samp[2];
		int u = samp[3];
		int tv = t-v;
		int st = s-t;
		int us = u-s;
		if(tv>=minsize && st>=minsize && us>=minsize){
			std::vector<double> vec_vt(t-v,0.0), vec_tu(u-t,0.0), vec_ts(s-t,0.0), vec_su(u-s,0.0);
			for(int i=v+1; i<=t; ++i)
				vec_vt[i-v-1] = Z[i];
			for(int i=t+1; i<=u; ++i)
				vec_tu[i-t-1] = Z[i];
			for(int i=t+1; i<=s; ++i)
				vec_ts[i-t-1] = Z[i];
			for(int i=s+1; i<=u; ++i)
				vec_su[i-s-1] = Z[i];
			kVec[j] = dist_ks(vec_vt,vec_tu) - dist_ks(vec_vt,vec_ts) - dist_ks(vec_ts,vec_su);
			++j;
		}
	}
	//Find 1-eps quantile estimate
	std::nth_element(kVec.begin(), kVec.begin() + (int)((1-eps)*R), kVec.end());
	return kVec[(int)((1-eps)*R)];
}

double dist_ks(std::vector<double>& X, std::vector<double>& Y){
	// Sorted X and Y, get their sizes.
	int n = X.size();
	int m = Y.size();
	// Maximum absolute difference between empirical
	// cdf values.
	double max_diff = 0;
	// Store the current difference between empirical cdf values.
	double current_diff = 0;
	// Index of element in X or Y currently under consideration.
	int x_counter = 0;
	int y_counter = 0;
	// Determine the maximum difference betwen the empirical
	// cdf's of X and Y.
	// Iterate over X and Y comparing the elements at index
	// x_counter and y_counter.
	while (x_counter < n && y_counter < m){
		if (X[x_counter] < Y[y_counter]){
			current_diff = current_diff + 1.0/n;
			++x_counter;
		} else if (X[x_counter] > Y[y_counter]){
			current_diff = current_diff - 1.0/m;
			++y_counter;		
		} else {
			// If values are equal then we need to increment both
			// the x_counter and y_counter.
			double s = X[x_counter];
			while (x_counter < n && X[x_counter] == s){
				current_diff = current_diff + 1.0/n;
				++x_counter;
			}
			while (y_counter < m && Y[y_counter] == s){
				current_diff = current_diff - 1.0/m;
				++y_counter;
			}
		}
		// Update the current value of max_diff.
		max_diff = (fabs(current_diff) > max_diff) ?
		fabs(current_diff) : max_diff;
	}
	// Since (x_counter < n && y_counter < m) is false then either
	// x_counter == n or y_counter == m. This while loop deals with the
	// first case. The following with the later.
	while (y_counter < m) {
		current_diff = current_diff - 1.0/m;
		++y_counter;
		max_diff = (fabs(current_diff) > max_diff) ?
		fabs(current_diff) : max_diff;
	}
	while (x_counter < n) {
		current_diff = current_diff + 1.0/n;
		++x_counter;
		max_diff = (fabs(current_diff) > max_diff) ?
		fabs(current_diff) : max_diff;
	}
	return max_diff;
}
