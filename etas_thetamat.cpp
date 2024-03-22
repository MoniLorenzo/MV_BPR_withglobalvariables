#include <Rcpp.h>
using namespace Rcpp;
using std::exp;

// [[Rcpp::export]]
NumericMatrix compute_log_etas_thetaMAT_cpp(bool Exp, NumericMatrix Linpred_mat, NumericMatrix THETA_mat) {
  int N = Linpred_mat.nrow();
  int Rym1 = Linpred_mat.ncol();
  
  //NumericMatrix thetamat;
  //if (THETA.nrow() == 1 && THETA.ncol() == 1) {
  //thetamat = NumericMatrix(N, Rym1, THETA(0, 0));
  //} else {
  //  thetamat = THETA;
  //}
  
  //NumericMatrix Nums_1_to_rym1(N, Rym1);
  //for (int i = 0; i < N; ++i) {
  //  for (int r = 0; r < Rym1; ++r) {
  //    Nums_1_to_rym1(i, r) = THETA(r) + Linpred_mat(i, r);
  //  }
  //}
  
  NumericMatrix Nums_1_to_ry(N, Rym1 + 1);
  for (int i = 0; i < N; ++i) {
    for (int r = 0; r < Rym1; ++r) {
      Nums_1_to_ry(i, r) = THETA_mat(i, r) + Linpred_mat(i, r);
    }
    Nums_1_to_ry(i, Rym1) = 0;
  }
  
  NumericMatrix etas_i(N, Rym1 + 1);
  for (int i = 0; i < N; ++i) {
    double max_val = Nums_1_to_ry(i, 0);
    for (int r = 1; r < Rym1 + 1; ++r) {
      if (Nums_1_to_ry(i, r) > max_val) {
        max_val = Nums_1_to_ry(i, r);
      }
    }
    double sum_exp = 0;
    for (int r = 0; r < Rym1 + 1; ++r) {
      etas_i(i, r) =  std::exp(Nums_1_to_ry(i, r) - max_val);
      sum_exp += etas_i(i, r);
    }
    for (int r = 0; r < Rym1 + 1; ++r) {
      etas_i(i, r) /= sum_exp;
    }
  }
  
  if (!Exp) {
    for (int i = 0; i < N; ++i) {
      for (int r = 0; r < Rym1 + 1; ++r) {
        etas_i(i, r) = log(etas_i(i, r));
      }
    }
  }
  return etas_i;
}