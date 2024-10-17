#include <Rcpp.h>
#include "MDP.h"

//#define DEBUG

using namespace Rcpp;

// Note: all are 0-based integer indices
// epsilon -1 means 0 for solved models and 1 for unsolved models

// [[Rcpp::export]]
NumericVector QV_cpp(
    NumericVector V,
    NumericMatrix P,
    NumericMatrix R,
    double GAMMA
    ) {
  
  int n = P.nrow();
  NumericVector V_prime(n);
  
  for (int s = 0; s < n; ++s) {
    V_prime[s] = sum(P(s, _) * (R(s, _) + GAMMA * V));
  }
 
  return(V_prime);
}