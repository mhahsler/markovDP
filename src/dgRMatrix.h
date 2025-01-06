#ifndef dGRMATRIX_H_
#define dGRMATRIX_H_

//  This is based on the dgCMatrix implementations at:
//  * https://gallery.rcpp.org/articles/sparse-matrix-class/
//  * https://codereview.stackexchange.com/questions/259594/rcpp-sparse-csc-matrix-class
//  * https://github.com/zdebruine/RcppSparse

#include <Rcpp.h>

namespace Rcpp {

class dgRMatrix {
public:
  IntegerVector i, p, Dim;
  NumericVector x;
  List Dimnames;
  
  // constructor
  dgRMatrix(S4 mat) {
    i = mat.slot("j");
    p = mat.slot("p");
    x = mat.slot("x");
    Dim = mat.slot("Dim");
    Dimnames = mat.slot("Dimnames");
  };
  
  int nrow() { return Dim[0]; };
  int ncol() { return Dim[1]; };
  int rows() { return Dim[0]; };
  int cols() { return Dim[1]; };
  
  double at(int row, int col) const {
    for (int j = p[row]; j < p[row + 1]; ++j) {
      if (i[j] == col) return x[j];
      else if (i[j] > col) break;
    }
    return 0.0;
  }
  double operator()(int row, int col) { return at(row, col); };
  
  NumericVector row(int row) const {
    NumericVector c(Dim[0], 0.0);
    for (int j = p[row]; j < p[row + 1]; ++j)
      c[i[j]] = x[j];
    return c;
  }
  
  NumericVector col(int col) const {
    NumericVector r(Dim[1], 0.0);
    for (int row = 0; row < Dim[1]; ++row) {
      for (int j = p[row]; j < p[row + 1]; ++j) {
        if (i[j] == col) r[row] = x[j];
        else if (i[j] > col) break;
      }
    }
    return r;
  }

  NumericMatrix dense() const {
    NumericMatrix res(Dim[0], Dim[1]);
    for (int j = 0; j < Dim[0]; ++j) {
      res.row(j) = row(j);
    }
    return res;
  }
    
};

// template <> dgRMatrix as(SEXP mat) { return dgRMatrix(mat); }
// 
// template <> SEXP wrap(const dgRMatrix& sm) {
//   S4 s(std::string("dgRMatrix"));
//   s.slot("i") = sm.i;
//   s.slot("p") = sm.p;
//   s.slot("x") = sm.x;
//   s.slot("Dim") = sm.Dim;
//   s.slot("Dimnames") = sm.Dimnames;
//   return s;
//}

}

#endif