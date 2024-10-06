#include <Rcpp.h>
#include "MDP.h"

//#define DEBUG

using namespace Rcpp;

// Note: all are 0-based integer indices
// epsilon -1 means 0 for solved models and 1 for unsolved models

// [[Rcpp::export]]
NumericVector QV_cpp(
    NumericVector U,
    NumericMatrix P,
    NumericMatrix R,
    double GAMMA
    ) {
  
  int n = P.nrow();
  NumericVector U_prime(n);
  
  for (int s = 0; s < n; ++s) {
    U_prime[s] = sum(P(s, _) * (R(s, _) + GAMMA * U));
  }
 
  return(U_prime);
}