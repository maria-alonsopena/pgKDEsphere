#include <RcppArmadillo.h>
#include <cmath>




using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat pgKDE_sph_vMF(arma::mat x, arma::mat datax,  double h, double cl, arma::vec mu, double nu){		
   
   int nx = x.n_rows;						// number of evaluation points
   int n = datax.n_rows;					// sample size
   arma:: mat kdex_aux(nx,n);				// auxiliary matrix: the ith row contains the KDE contribution for evaluation point xi
   arma::vec kde(nx);						// aux vector: the KDE (without constants) for each evaluation point
   for(int i=0; i<nx; ++i){
	   
	   arma::mat xxi = x.row(i);
	   arma::mat aux = exp(-(1-datax * xxi.t())/(h*h)) * exp(nu *  xxi * mu)/ exp(nu * datax * mu);
	   kdex_aux.row(i) = aux.t();
	   kde(i) = sum( kdex_aux.row(i) );

	}
   
    return kde/(cl*n);

}




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat KDE_sph_vMF(arma::mat x, arma::mat datax,  double h, double cl){		
   
   int nx = x.n_rows;						// number of evaluation points
   int n = datax.n_rows;					// sample size
   arma:: mat kdex_aux(nx,n);				// auxiliary matrix: the ith row contains the KDE contribution for evaluation point xi
   arma::vec kde(nx);						// aux vector: the KDE (without constants) for each evaluation point
   for(int i=0; i<nx; ++i){
	   
	   arma::mat xxi = x.row(i);
	   arma::mat aux = exp(-(1-datax * xxi.t())/(h*h));
	   kdex_aux.row(i) = aux.t();
	   kde(i) = sum( kdex_aux.row(i) );

	}
   
    return kde/(cl*n);

}


arma::colvec Arma_rowSums(const arma::mat& x) {
  return arma::sum(x, 1);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec KDE(arma::mat x, arma::mat datax, double h, double cl){		
   
   int nx = x.n_rows;
   int n = datax.n_rows;
   double hh = h*h;
   arma::mat vecX;
   arma::vec kde(nx);
   double sumx;
   
   for(int i=0; i<nx; ++i){
	   
	   vecX = exp(-(1-x.row(i) * datax.t())/hh);
	   sumx = sum(Arma_rowSums(vecX));
	   kde(i) = sumx;
	
   }

	
    return kde/(cl * n);

}
