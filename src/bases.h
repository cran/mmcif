#ifndef SPLINES_H
#define SPLINES_H

#include "arma-wrap.h"
#include <limits> // quiet_NaN
#include <stdexcept> // invalid_argument
#include <memory>
#include "wmem.h"
#include "lp-mmcif.h"
#include <algorithm>

namespace bases {

constexpr unsigned default_order{4};
constexpr int default_ders{0};
constexpr bool default_intercept{false},
                 default_use_log{false};

using namespace arma;

/// base class for basis expansions
class basisMixin {
public:
  /// the required working memory
  virtual size_t n_wmem() const = 0;

  /// the number of basis functions
  virtual unsigned n_basis() const = 0;

  /**
   * fills a vector with the (possibly derivatives) of the basis expansions
   * evaluated at x */
  void operator()
    (vec &out, double *wk_mem, double const x,
     const int ders = default_ders) const {
    (*this)(out.memptr(), wk_mem, x, ders);
  }
  /// returns an allocated vector
  vec operator()
    (double const x, double * wk_mem,
     int const ders = default_ders) const {
    vec out(n_basis());
    (*this)(out.begin(), wk_mem, x, ders);
    return out;
  }
  /// same as the other operator() calls but filling the out
  virtual void operator()
    (double *out, double *wk_mem, double const x,
     int const ders = default_ders) const = 0;

  mat basis
    (const vec &x, double *wk_mem, const int ders = default_ders,
     const double centre = std::numeric_limits<double>::quiet_NaN())
    const  {
    unsigned const n_basis_v(n_basis()),
    n_x    (x.n_elem);
    rowvec centering =
      (std::isnan(centre) || ders > 0 ?
      zeros(n_basis_v) : operator()(centre, wk_mem, 0)).t();

    mat out(n_x, n_basis_v);
    vec wrk(n_basis_v);
    for (unsigned i = 0; i < n_x; i++){
      (*this)(wrk, wk_mem, x[i], ders);
      out.row(i) = wrk.t() - centering;
    }

    return out;
  }

  virtual ~basisMixin() = default;

  virtual std::unique_ptr<basisMixin> clone() const = 0;

  virtual void set_lower_limit(double const x){
    if(use_log)
      lower_limit = log(x);
    else
      lower_limit = x;
  }

  basisMixin(bool const use_log = default_use_log):
    use_log{use_log} {
      set_lower_limit(use_log ? std::numeric_limits<double>::epsilon() : 0);
    }

protected:
  /// should the log of the input be used
  bool use_log;
  /// lower limit of integrals. This is log transformed if use_log is true
  double lower_limit;
};

class SplineBasis : public basisMixin {
    void comp_basis(double const x, double *out,
                  double * wk_mem, unsigned const ders) const {
    // Cox-de Boor recursion formula
    //    https://en.wikipedia.org/wiki/De_Boor%27s_algorithm

    // find the knot such that knot[i] <= x < knot[i + 1]
    double const * it_inter{std::upper_bound(knots.begin(), knots.end(), x)};

    // deal with the matching boundaries
    bool const is_at_boundary{it_inter == knots.end() && *std::prev(it_inter) == x};
    while((it_inter == knots.end() && *std::prev(it_inter) == x) ||
          (it_inter != knots.end() && *it_inter == x &&
           it_inter != knots.begin() && *it_inter == *std::prev(it_inter)))
      --it_inter;

    std::fill(out, out + ncoef, 0);
    if(it_inter == knots.begin() || it_inter == knots.end() || ordm1 < ders)
      return;

    --it_inter;
    unsigned const shift = std::distance(knots.begin(), it_inter);
    double * const D{wk_mem};

    // set the initial one
    std::fill(D, D + order, 0);
    D[ordm1] = !is_at_boundary || ordm1 > ders;

    bool const is_interior{ordm1 <= shift && shift + ordm1 < knots.size()};
    if(is_interior){
      for(unsigned r = 1; r <= ordm1 - ders; ++r){
        unsigned const j_start{order - r};
        double const *k1{knots.begin() + shift + j_start - ordm1},
                     *k2{k1 + r};

        for(unsigned j = j_start; j < order; ++j, ++k1, ++k2){
          double const w_new{*k1 == *k2 ? 0 : (x - *k1) / (*k2 - *k1)};

          // update the previous
          D[j - 1] += (1 - w_new) * D[j];
          // update this one
          D[j] *= w_new;
        }
      }

      // handle the derivatives
      for(unsigned r = ordm1 - ders + 1; r <= ordm1; ++r){
        unsigned const j_start{order - r};
        double const *k1{knots.begin() + shift + j_start - ordm1},
                     *k2{k1 + r};

        for(unsigned j = j_start; j < order; ++j, ++k1, ++k2){
          double const w_new{*k1 == *k2 ? 0 : r / (*k2 - *k1)};

          // update the previous
          D[j - 1] -= w_new * D[j];
          // update this one
          D[j] *= w_new;
        }
      }

    } else {
      unsigned const j_min{ordm1 > shift ? ordm1 - shift : 1};

      for(unsigned r = 1; r <= ordm1 - ders; ++r){
        unsigned const j_start{std::max(order - r, j_min)},
                       j_max  {std::min(order, knots.size() - shift + ordm1 - r)};
        double const *k1{knots.begin() + shift + j_start - ordm1},
                     *k2{k1 + r};

        for(unsigned j = j_start; j < j_max; ++j, ++k1, ++k2){
          double const w_new{*k1 == *k2 ? 0 : (x - *k1) / (*k2 - *k1)};

          // update the previous
          D[j - 1] += (1 - w_new) * D[j];
          // update this one
          D[j] *= w_new;
        }
      }

      // handle the derivatives
      for(unsigned r = ordm1 - ders + 1; r <= ordm1; ++r){
        unsigned const j_start{std::max(order - r, j_min)},
                       j_max  {std::min(order, knots.size() - shift + ordm1 - r)};
        double const *k1{knots.begin() + shift + j_start - ordm1},
                     *k2{k1 + r};

        for(unsigned j = j_start; j < j_max; ++j, ++k1, ++k2){
          double const w_new{*k1 == *k2 ? 0 : r / (*k2 - *k1)};

          // update the previous
          D[j - 1] -= w_new * D[j];
          // update this one
          D[j] *= w_new;
        }
      }
    }

    if(shift < ordm1){
      unsigned const shift_D{ordm1 - shift};
      std::copy(D + shift_D, D + order, out);
    } else {
      unsigned const shift_out{shift - ordm1},
                          n_cp{std::min(ncoef - shift_out, order)};
      std::copy(D, D + n_cp, out + shift_out);
    }
  }

public:
  unsigned const order = default_order, /* order of the spline */
                     ordm1 = order - 1;     /* order - 1 (3 for cubic splines) */
  vec const knots;	               /* knot vector */
  unsigned const nknots = knots.n_elem, /* number of knots
                                           except for the boundary case */
                      ncoef =               /* number of coefficients */
                         nknots > order ? nknots - order : 0L;

  SplineBasis(const vec &knots, const unsigned order = default_order,
              bool const use_log = default_use_log,
              bool const with_integral = true);

  SplineBasis(SplineBasis const &other):
    SplineBasis(other.knots, other.order, other.use_log,
                static_cast<bool>(other.integral_basis)) { }

  unsigned n_basis() const override {
    return ncoef;
  }

  size_t n_wmem() const override {
    return n_wmem_v;
  }

  using basisMixin::set_lower_limit;
  using basisMixin::operator();

  /** the function ignores the use_log as this class is usually not used on its
   * own
   */
  void operator()
    (double *out, double *wk_mem, double const x,
     const int ders = default_ders)
    const override {
    if(ders >= 0){
      comp_basis(x, out, wk_mem, ders);
      return;
    }

    if(ders < -1)
      throw std::runtime_error("not implemented for ders < -1");
    // use formulas from Monotone Regression Splines in Action
    //    https://doi.org/10.1214/ss/1177012761
    // TODO: can be implemented smarter...
    double * const basis_mem = wk_mem;
    wk_mem += integral_basis->n_basis();

    double dorder{static_cast<double>(order)};

    // computes the indefinte integral at the upper or lower limit. The
    // function must first be called at the upper limit
    auto add_int = [&](double lim, bool const is_upper){
      // we may integrate up to the limit but no further
      lim = std::min(knots.back(), lim);

      // evaluate the basis which is one order greater
      (*integral_basis)(basis_mem, wk_mem, lim, ders + 1);

      // find the index j such that knots[j] <= lim < knots[j + 1]
      // use -1 if x < knots[0]
      auto const idx_knot_start =
        ([&]{
          auto knot_j = std::upper_bound(knots.begin(), knots.end(), lim);
          return static_cast<unsigned>
            (std::distance(knots.begin(), knot_j) - 1);
        })();

      // x is too small for these basis function to be active
      unsigned const idx_no_support
        {std::min<unsigned>(idx_knot_start + 1, ncoef)};
      if(is_upper)
        std::fill(out + idx_no_support, out + ncoef, 0);
      // x is large enough that we have integrated over the full support of
      // these basis functions
      unsigned i{};
      unsigned const is_capped =
        idx_knot_start + 1 >= order ? idx_knot_start + 1 - order : 0;
      for(; i < is_capped; ++i)
        if(is_upper)
          out[i]  = (knots[i + order] - knots[i]) / dorder;
        else
          out[i] -= (knots[i + order] - knots[i]) / dorder;

      // the residual is somewhere in between
      for(; i < idx_no_support; ++i){
        double su{};
        for(unsigned j = i; j < idx_no_support; ++j)
          // TODO: redundant computations
          su += basis_mem[j];
        if(is_upper)
          out[i]  = su * (knots[i + order] - knots[i]) / dorder;
        else
          out[i] -= su * (knots[i + order] - knots[i]) / dorder;
      }
    };

    add_int(x, true); // the upper limit
    if(lower_limit > knots[0])
      add_int(lower_limit, false); // the lower limit
  }

  virtual ~SplineBasis() = default;

  std::unique_ptr<basisMixin> clone() const override {
    return std::make_unique<SplineBasis>(*this);
  }

private:
  std::unique_ptr<SplineBasis> integral_basis;
  size_t const n_wmem_v
    {
    integral_basis
      ? integral_basis->n_wmem() + integral_basis->n_basis()
      : ordm1
    };
};

class bs final : public SplineBasis {
  void do_eval
    (double *out, double *wk_mem, double const x,
     const int ders = default_ders) const {
    double * const my_wk_mem{wk_mem};
    wk_mem += std::max(SplineBasis::n_basis(), bs::n_basis());

    if(ders < 0){
      // handle integration. First we handle the interior part. Then we address
      // the extrapolation if needed
      if(intercept)
        SplineBasis::operator()(out, wk_mem, x, ders);
      else {
        SplineBasis::operator()(my_wk_mem, wk_mem, x, ders);
        for(unsigned i = 1; i < SplineBasis::n_basis(); ++i)
          out[i - 1L] = my_wk_mem[i];
      }

      auto handle_outside = [&](double const x, double const sign){
        if(x >= boundary_knots[0] && x <= boundary_knots[1])
          return;

        double const k_pivot =
          x < boundary_knots[0]
            ? 0.75 * boundary_knots[0] + 0.25 * knots[order]
            : 0.75 * boundary_knots[1] + 0.25 * knots[knots.n_elem - order - 2],
                       delta = x - k_pivot,
                 delta_bound =
          x < boundary_knots[0]
            ? boundary_knots[0] - k_pivot
            : boundary_knots[1] - k_pivot;

        auto add_term = [&](int const d, double const f = 1){
          bs::operator()(my_wk_mem, wk_mem, k_pivot, d);
          for(unsigned i = 0; i < bs::n_basis(); ++i)
            out[i] += sign * f * my_wk_mem[i];
        };

        double m1{1}, m2{1}, denom{1};
        for(unsigned i = 1; i <= 4; ++i){
          m1 *= delta;
          m2 *= delta_bound;
          denom *= i;
          add_term(i - 1, (m1 - m2) / denom);
        }
      };

      handle_outside(x          , 1);
      handle_outside(lower_limit, -1);
      return;
    }

    if(x < boundary_knots[0] || x > boundary_knots[1]) {
      double const k_pivot =
        x < boundary_knots[0]
          ? 0.75 * boundary_knots[0] + 0.25 * knots[order]
          : 0.75 * boundary_knots[1] + 0.25 * knots[knots.n_elem - order - 2],
                     delta = x - k_pivot;

      auto add_term = [&](int const d, double const f = 1){
        do_eval(my_wk_mem, wk_mem, k_pivot, d);
        for(unsigned i = 0; i < bs::n_basis(); ++i)
          out[i] += f * my_wk_mem[i];
      };

      std::fill(out, out + bs::n_basis(), 0);

      add_term(ders);
      double m1{1};
      for(unsigned i = ders + 1; i < order; ++i){
        m1 *= delta / (i - ders);
        add_term(i, m1);
      }

      return;
    }

    if(intercept)
      SplineBasis::operator()(out, wk_mem, x, ders);
    else {
      SplineBasis::operator()(my_wk_mem, wk_mem, x, ders);
      for(unsigned i = 1; i < SplineBasis::n_basis(); ++i)
        out[i - 1L] = my_wk_mem[i];
    }
  }

public:
  double const boundary_knots[2];
  bool const intercept;
  unsigned const df;
  size_t n_wmem_v
    {2 * std::max(SplineBasis::n_basis(), bs::n_basis()) +
      SplineBasis::n_wmem()};

  size_t n_wmem() const {
    return n_wmem_v;
  }

  bs(const vec &bk, const vec &ik,
     const bool inter = default_intercept,
     const unsigned ord = default_order,
     bool const use_log = default_use_log);

  unsigned n_basis() const {
    return SplineBasis::n_basis() - (!intercept);
  }

  std::unique_ptr<basisMixin> clone() const {
    return std::make_unique<bs>(*this);
  }

  using SplineBasis::set_lower_limit;
  using SplineBasis::operator();
  void operator()
      (double *out, double *wk_mem, double const x,
       const int ders = default_ders) const {
      if(!use_log){
        do_eval(out, wk_mem, x, ders);
        return;
      }

      do_eval(out, wk_mem, std::log(x), ders);
      switch(ders){
      case 0:
        break;
      case 1:
        for(unsigned i = 0; i < n_basis(); ++i)
          // TODO: we can use that some are zero
          out[i] /= x;
        break;
      default:
        throw std::runtime_error
        ("not implemented with use_log and ders " + std::to_string(ders));
    }
  }
};

class ns final : public basisMixin {
  SplineBasis s_basis;

  void do_eval
    (double *out, double *wk_mem, double const x,
     int const ders = default_ders) const {
    if(ders < 0){
      if(ders < -1)
        throw std::runtime_error("integration not implemented for order 2 or higher");

      // let K0 be the smallest knot. Then the spline does the right thing as
      // long as K0 <= lb <= KMax and and K0 <= ub <= KMax. Otherwise, we have
      // to make adjustments

      // handle the integration between K0 and KMax
      {
        double * const lhs = wk_mem;
        wk_mem += q_matrix.n_rows;
        double * const b = wk_mem;
        wk_mem += s_basis.n_basis();
        s_basis(b, wk_mem, x, ders);

        std::fill(lhs, lhs + q_matrix.n_rows, 0);
        lp_mmcif::mat_vec
          (lhs, q_matrix.begin(), b + (!intercept), q_matrix.n_rows,
           q_matrix.n_cols);

        std::copy(lhs + 2, lhs + q_matrix.n_rows, out);
      }

      // handle the areas outside of the knots
      auto handle_outside = [&](double const x, double const sign){
        if(x < boundary_knots[0]){
          double const b{boundary_knots[0]};
          for(unsigned i = 0; i < ns::n_basis(); ++i){
            double const new_term{
              tl1[i] * x * (x / 2 - b) + x * tl0[i]  -
                tl1[i] * b * (b / 2 - b) - b * tl0[i]};
            out[i] += sign * new_term;
          }

        } else if(x > boundary_knots[1]){
          double const b{boundary_knots[1]};
          for(unsigned i = 0; i < ns::n_basis(); ++i){
            double const new_term
              {tr1[i] * x * (x / 2 - b) + x * tr0[i] -
                tr1[i] * b * (b / 2 - b) - b * tr0[i]};
            out[i] += sign * new_term;
          }

        }
      };

      handle_outside(x          , 1);
      handle_outside(lower_limit, -1);
      return;
    }

    // ders >= 0
    if(x < boundary_knots[0]){
      if(ders==0){
        for(unsigned i = 0; i < ns::n_basis(); ++i){
          out[i] = tl1[i];
          out[i] *= x - boundary_knots[0];
          out[i] += tl0[i];
        }

      } else if (ders == 1)
        std::copy(tl1.begin(), tl1.end(), out);
      else
        std::fill(out, out + ns::n_basis(), 0);

      return;

    } else if (x > boundary_knots[1]) {
      if (ders==0){
        for(unsigned i = 0; i < ns::n_basis(); ++i){
          out[i] = tr1[i];
          out[i] *= x - boundary_knots[1];
          out[i] += tr0[i];
        }

      } else if (ders==1)
        std::copy(tr1.begin(), tr1.end(), out);
      else
        std::fill(out, out + ns::n_basis(), 0);

      return;
    }

    double * const lhs = wk_mem;
    wk_mem += q_matrix.n_rows;
    double * const b = wk_mem;
    wk_mem += s_basis.n_basis();
    s_basis(b, wk_mem, x, ders);

    std::fill(lhs, lhs + q_matrix.n_rows, 0);
    lp_mmcif::mat_vec
      (lhs, q_matrix.begin(), b + (!intercept), q_matrix.n_rows,
       q_matrix.n_cols);

    std::copy(lhs + 2, lhs + q_matrix.n_rows, out);
  }

  vec trans(const vec &x) const {
    // TODO: very inefficient
    vec out = q_matrix * (intercept ? x : x(span(1, x.n_elem - 1)));
    return out(span(2, out.size() - 1));
  }

public:
  double const boundary_knots[2];
  bool const intercept;
  mat const q_matrix = ([&](){
    // calculate the Q matrix
    arma::vec tmp{boundary_knots[0], boundary_knots[1]};
    mat const_basis = s_basis.basis
      (tmp, wmem::mem_stack().get(s_basis.n_wmem()), 2);
    if (!intercept)
      const_basis = const_basis.cols(1, const_basis.n_cols - 1);
    mat qd, rd;
    if(!qr(qd, rd, const_basis.t()))
      throw std::invalid_argument("ns: QR decomposition failed");
    inplace_trans(qd);
    return qd;
  })();
  vec const tl0, tl1, tr0, tr1;

  ns(const vec &boundary_knots, const vec &interior_knots,
     const bool intercept = default_intercept,
     const unsigned order = default_order,
     const bool use_log = default_use_log);

  size_t n_wmem() const {
    return s_basis.n_wmem() + q_matrix.n_rows +
      s_basis.n_basis() + n_basis();
  }

  unsigned n_basis() const {
    return q_matrix.n_rows - 2;
  }

  std::unique_ptr<basisMixin> clone() const {
    return std::make_unique<ns>(*this);
  }

  void set_lower_limit(double const x){
    basisMixin::set_lower_limit(x);
    s_basis.set_lower_limit(x);
  }

  using basisMixin::operator();
  void operator()
    (double *out, double *wk_mem, double const x,
     int const ders = default_ders) const {
    if(!use_log){
      do_eval(out, wk_mem, x, ders);
      return;
    }

    do_eval(out, wk_mem, std::log(x), ders);
    switch(ders){
    case 0:
      break;
    case 1:
      for(unsigned i = 0; i < n_basis(); ++i)
        // TODO: we can use that some are zero
        out[i] /= x;
      break;
    default:
      throw std::runtime_error
      ("not implemented with use_log and ders " + std::to_string(ders));
    }
  }
}; // class ns

class iSpline final : public basisMixin {
public:
  bool const intercept;
  unsigned const order;
  bs bspline; // TODO: can be a SplineBasis

public:
  iSpline(const vec &boundary_knots, const vec &interior_knots,
          const bool intercept = default_intercept,
          const unsigned order = default_order);

  size_t n_wmem() const {
    return bspline.n_wmem() + bspline.n_basis();
  }

  unsigned n_basis() const {
    return bspline.n_basis() - (!intercept);
  }

  std::unique_ptr<basisMixin> clone() const {
    return std::make_unique<iSpline>(*this);
  }

  void set_lower_limit(double const x){
    basisMixin::set_lower_limit(x);
    bspline.set_lower_limit(x);
  }

  using basisMixin::operator();
  void operator()
    (double *out, double *wk_mem, double const x,
     int const ders = default_ders) const  {
    double * const b{wk_mem};
    unsigned const n_b{bspline.n_basis()};
    wk_mem += n_b;

    if(x < 0){
      std::fill(out, out + n_basis(), 0);
      return;

    }
    else if(x <= 1){
      bspline(b, wk_mem, x, ders);
      unsigned const js = bspline.knots.size() - 2 > 0 ?
        static_cast<unsigned>(std::lower_bound(
          bspline.knots.begin(),
          /* TODO: should this not be end and not -1? */
          bspline.knots.end() - 1L, x) -
          bspline.knots.begin()) :
        order + 1;
      for(unsigned j = n_b; j-- >0;)
        if(j > js)
          b[j] = 0.0;
        else if(j != n_b - 1)
          b[j] += b[j+1];
      if(ders==0)
        for(unsigned j = n_b - 1; j-- > 0;)
          if(j + order + 1 < js)
            b[j] = 1.0;

      if(intercept)
        std::copy(b, b + bspline.n_basis(), out);
      else
        std::copy(b + 1, b + bspline.n_basis(), out);
      return;

    }
    else if(ders > 0)
      std::fill(out, out + n_basis(), 0);
    else
      std::fill(out, out + n_basis(), 1);
  }
}; // class iSpline

class mSpline final : public basisMixin {
public:
  bs bspline; // TODO: can be a SplineBasis
  bool const intercept;

public:
  mSpline(const vec &boundary_knots, const vec &interior_knots,
          const bool intercept = default_intercept,
          const unsigned order = default_order);

  size_t n_wmem() const {
    return bspline.n_wmem() + bspline.n_basis();
  }

  unsigned n_basis() const {
    return bspline.n_basis() - (!intercept);
  }

  std::unique_ptr<basisMixin> clone() const {
    return std::make_unique<mSpline>(*this);
  }

  void set_lower_limit(double const x){
    basisMixin::set_lower_limit(x);
    bspline.set_lower_limit(x);
  }

  using basisMixin::operator();
  void operator()
    (double *out, double *wk_mem, double const x,
     int const ders = default_ders) const {
    double * const wrk{wk_mem};
    wk_mem += bspline.n_basis();

    bspline(wrk, wk_mem, x, ders);
    for (unsigned j = 0; j < bspline.n_basis(); ++j) {
      double denom = bspline.knots(j + bspline.order) - bspline.knots(j);
      wrk[j] *= denom > 0.0 ? bspline.order / denom : 0.0;
    }

    if(intercept)
      std::copy(wrk, wrk + bspline.n_basis(), out);
    else
      std::copy(wrk + 1, wrk + bspline.n_basis(), out);
  }
}; // class mSpline

class orth_poly final : public basisMixin {
  // coefficients for the orthogonal polynomial
  vec alpha,
      norm2,
      sqrt_norm2{arma::sqrt(norm2)};
  // flags for whether a raw polynomial is used and whether there is a intercept
  bool raw,
       intercept;
  // the number of basis function plus the possible intercept
  unsigned n_basis_v;
  // the matrix to map from the raw polynomial to the orthogonal polynomial
  // see https://stats.stackexchange.com/a/472289/81865
  std::vector<double> orth_map;

  // evaluates the polynomial with raw == TRUE
  static void eval_raw
    (double *out, double const x, bool const inter, int const ders,
     unsigned const degree, double const lb) {
    unsigned const dim{degree + inter};

    if(ders == 0){
      double val{inter ? 1 : x};
      for(unsigned c = 0; c < dim; ++c, val *= x)
        out[c] = val;
      return;

    } else if(ders < 0){
      // compute the starting value
      double val_upper{x},
             val_lower{lb};
      unsigned const uders{static_cast<unsigned>(-ders)};
      for(unsigned i = 2; i <= uders; ++i){
        val_upper *= x  / static_cast<double>(i);
        val_lower *= lb / static_cast<double>(i);
      }

      if(!inter){
        val_upper *= x  / (uders + 1);
        val_lower *= lb / (uders + 1);
      }

      for(unsigned c = 0; c < dim; c++){
        out[c] = val_upper - val_lower;
        val_upper *= x  / static_cast<double>(c + uders + 1 + !inter);
        val_lower *= lb / static_cast<double>(c + uders + 1 + !inter);
        if(c + 1 + !inter >= uders){
          val_upper *= c + 1. + !inter;
          val_lower *= c + 1. + !inter;
        }
      }

    } else { // ders > 0
      unsigned const uders{static_cast<unsigned>(ders)};

      if(inter){
        std::fill(out, out + uders, 0);
        double val{1};
        for(unsigned c = uders; c < dim; c++){
          unsigned mult{c};
          for(unsigned cc = c; --cc > c - uders;)
            mult *= cc;
          out[c] = mult * val;
          val *= x;
        }

      } else {
        std::fill(out, out + uders - 1, 0);
        double val{1};
        for(unsigned c = uders - 1; c < dim; c++){
          unsigned mult{c + 1};
          for(unsigned cc = c + 1; --cc > c - uders + 1;)
            mult *= cc;
          out[c] = mult * val;
          val *= x;
        }
      }
    }
  }

  void do_eval(double *out, double *wk_mem, double const x,
               int const ders)
    const {
    if(raw){
      eval_raw(out, x, intercept, ders, n_basis_v - intercept, lower_limit);
      return;
    }

    if(ders == 0){
      out[0] = 1.;
      double old{1};

      if(alpha.n_elem > 0L){
        out[0 + intercept] = x - alpha[0];
        for(unsigned c = 1; c < alpha.n_elem; c++){
          out[c + intercept] =
            (x - alpha[c]) * out[c - 1 + intercept] - norm2[c + 1] /
            norm2[c] * old;
          old = out[c - 1 + intercept];
        }
      }

      for(unsigned j = 1; j < alpha.n_elem + 1; ++j)
        out[j - 1 + intercept] /= sqrt_norm2.at(j + 1);

      return;
    }

    // compute the raw polynomial and multiply on the matrix
    // TODO: can likely be done in more stable way?
    eval_raw(wk_mem, x, true, ders, n_basis_v - intercept, lower_limit);

    std::fill(out, out + n_basis_v, 0);
    // handle the intercept term
    auto g = orth_map.begin() + !intercept;
    if(intercept)
      out[0] = *g++ * wk_mem[0];

    // handle the other terms
    for(unsigned j = 1; j < alpha.size() + 1; ++j)
      for(unsigned i = 0; i <= j; ++i)
        out[j - !intercept] += wk_mem[i] * *g++;
  }

public:
  // constructor for raw == true
  orth_poly
    (unsigned const degree, bool const intercept = default_intercept,
     bool const use_log = default_use_log);

  // constructor for raw == false
  orth_poly(vec const &alpha, vec const &norm2,
            bool const intercept = default_intercept,
            bool const use_log = default_use_log);

  std::unique_ptr<basisMixin> clone() const {
    return std::make_unique<orth_poly>(*this);
  }

  size_t n_wmem() const {
    return intercept ? n_basis_v : n_basis_v + 1;
  }

  /**
   * behaves like predict(<poly object>, newdata) except there might be an
   * intercept */
  using basisMixin::set_lower_limit;
  using basisMixin::operator();
  void operator()(double *out, double *wk_mem, double const x,
                  int const ders = default_ders)
    const {
    if(!use_log){
      do_eval(out, wk_mem, x, ders);
      return;
    }

    do_eval(out, wk_mem, std::log(x), ders);
    switch(ders){
    case 0:
      break;
    case 1:
      for(unsigned i = 0; i < n_basis_v; ++i)
        out[i] /= x;
      break;
    default:
      throw std::runtime_error
          ("not implemented with use_log and ders " + std::to_string(ders));
    }
  }

  /**
   * behaves like poly(x, degree). The orthogonal polynomial is returned by
   * reference.
   */
  static orth_poly poly_basis(vec, unsigned const, mat&);

  unsigned n_basis() const {
    return n_basis_v;
  }
}; // orth_poly

using bases_vector = std::vector<std::unique_ptr<basisMixin> >;

/// simple function clone bases
inline bases_vector clone_bases(const bases_vector &bases){
  std::vector<std::unique_ptr<basisMixin> > out;
  out.reserve(bases.size());
  for(auto &x : bases)
    out.emplace_back(x->clone());
  return out;
}

} // namespace bases

#endif // SPLINES_H
