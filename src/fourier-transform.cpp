#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::cx_cube fsd_rcpp_fourier(const arma::cx_cube & A, const arma::mat & lags, const arma::mat & freq) {

  int d1 = A.n_rows;
  int d2 = A.n_cols;

  int freqlength = freq.n_cols;
  int lagslength = lags.n_cols;

  arma::mat lagst = lags.t();

  arma::cx_cube B(d1, d2, freqlength, fill::zeros);
  arma::cx_vec expweight(lagslength);

  arma::cx_double ii(0,1);

  for (int j0 = 0; j0 < freqlength; j0++) {
    expweight = arma::exp( ii * (lagst * freq.col(j0)) );
    for (int k = 0; k < lagslength; k++) {
      B.slice(j0) = B.slice(j0) + A.slice(k) * expweight(k);
    }
  }

  return B;
}

// [[Rcpp::export]]
arma::cx_cube fsd_rcpp_fourier_inverse(const arma::cx_cube & B, const arma::mat & freq, const arma::mat & lags) {
  int freqlength = freq.n_cols;
  return fsd_rcpp_fourier(B, freq, -lags) / freqlength;
}
