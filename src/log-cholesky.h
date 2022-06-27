#ifndef LOG_CHOLESKY_H
#define LOG_CHOLESKY_H
#include "arma-wrap.h"
#include "simple-mem-stack.h"

namespace log_chol {
struct pd_mat {
  /**
   * return the required memory to get the original matrix from a log Cholesky
   * decomposition */
  static size_t n_wmem(size_t const dim){
    return dim * dim;
  }

  /** Computes L^TL where L is a upper triangular matrix. The argument is a
    * a vector with the non-zero elements in column major order. The diagonal
    * entries are on the log scale. The last element is working memory */
  static void get
    (double const *theta, size_t const dim, double * res, double *wk_mem){
    arma::mat L(wk_mem, dim, dim, false, true);
    L.zeros();

    for(size_t j = 0; j < dim; ++j){
      for(size_t i = 0; i < j; ++i)
        L.at(i, j) = *theta++;
      L.at(j, j) = std::exp(*theta++);
    }

    arma::mat res_mat(res, dim, dim, false, true);
    res_mat = L.t() * L;
  }

  /// same as the above but perform the allocation of working memory
  static void get
    (double const *theta, size_t const dim, double * res,
     ghqCpp::simple_mem_stack<double> &mem_stack){
    get(theta, dim, res, mem_stack.get(n_wmem(dim)));
  }
};

struct dpd_mat {
  /**
   * return the required memory to get the derivative as part of the chain
   * rule  */
  static size_t n_wmem(size_t const dim){
    return 3 * dim * dim;
  }

  /**
   * computes the derivative w.r.t. theta as part of the chain rule. That is,
   * the derivatives of f(L^top L) where d/dX f(X) evaluated at X = L^top L
   * is supplied. Only the derivatives w.r.t. the upper triangular of X need
   * to be valid.
   */
  static void get
    (double const *theta, size_t const dim, double * __restrict__ res,
     double const * derivs, double * __restrict__ wk_mem){
    arma::mat L(wk_mem, dim, dim, false, true);
    L.zeros();

    for(size_t j = 0; j < dim; ++j){
      for(size_t i = 0; i < j; ++i)
        L.at(i, j) = *theta++;
      L.at(j, j) = std::exp(*theta++);
    }

    // TODO: can be done more efficient
    arma::mat dum(const_cast<double*>(derivs), dim, dim, false, true),
                D(L.end(), dim, dim, false, true);
    D = arma::symmatu(dum);
    arma::mat inter(D.end(), dim, dim, false, true);
    inter = L * D;

    double * __restrict__ r = res;
    for(size_t j = 0; j < dim; ++j){
      for(size_t i = 0; i < j; ++i)
        *r++ += 2 * inter.at(i, j);
      *r++ += 2 * inter.at(j, j) * L.at(j, j);
    }
  }

  /// same as the above but perform the allocation of working memory
  static void get
    (double const *theta, size_t const dim, double * __restrict__ res,
     double const * derivs, ghqCpp::simple_mem_stack<double> &mem_stack){
    get(theta, dim, res, derivs, mem_stack.get(n_wmem(dim)));
  }
};
} // namespace log_chol

#endif
