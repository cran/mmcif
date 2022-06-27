#include <testthat.h>
#include "lp-mmcif.h"
#include <vector>

context("the lp_mmcif functions work") {
  test_that("general_lsolver works") {
    /*
     n <- 5
     nrhs <- 3
     set.seed(777)
     dput(X <- rnorm(n * n) |> matrix(n) |> round(3))
     dput(z <- rnorm(nrhs * n) |> round(3) |> matrix(n))
     dput(solve(X, z))
     dput(solve(t(X), z))
     */
    constexpr int n{5}, nrhs{3};
    constexpr double A[]{0.49, -0.399, 0.511, -0.399, 1.639, 0.621, 0.203, 1.109, -0.206, -0.379, -0.304, 0.054, -1.881, -0.034, 2.311, 0.972, 0.965, -0.544, 0.671, 0.501, -2.031, 0.228, -0.783, 1.274, 1.438},
                     z[]{.427, -1.744, -0.025, -1.488, -0.542, 0.661, 0.858, 1.252, -1.215, 1.523, -1.236, 1.666, -0.504, -0.578, 0.272},
                 truth[]{0.963248602818228, -2.453962925843, -0.594910953539176, -0.63648881984902, -0.943738731947272, -1.32338541800779, 5.86620682942923, 2.65522887230908, -1.09447298779662, 0.227689938420218, -3.62671282498382, 7.42283417643874, 3.8964949339858, -1.69633089489341, 0.608135228733778},
               truth_T[]{1.55719468317323, -5.50156649129805, -1.0589119957648, 2.9505451976936, -0.495899244092278, -2.49543779925809, 4.06633426334945, 1.00562643605077, -3.88551368459928, 0.879827226196364, -2.82629641610645, 6.16989418662416, 1.08795551138814, -4.81865155295422, 0.0805892546718268};

    std::vector<double> mem(lp_mmcif::general_lsolver::n_wk_mem(n));
    lp_mmcif::general_lsolver solver(n, A, mem.data());

    double res[n * nrhs];
    solver(z, res, nrhs, false);
    for(int i = 0; i < n * nrhs; ++i)
      expect_true(std::abs(res[i] - truth[i]) < std::abs(truth[i]) * 1e-8);

    solver(z, res, nrhs, true);
    for(int i = 0; i < n * nrhs; ++i)
      expect_true(std::abs(res[i] - truth_T[i]) < std::abs(truth_T[i]) * 1e-8);
  }

  test_that("backprop_cond_mean works as expected"){
    /*
     n <- 9L
     k1 <- 2L
     k2 <- 3L
     l1 <- 5L
     l2 <- 8L
     set.seed(22)
     dput(Sig <- rWishart(1, n, diag(n)) |> drop() |> round(2))
     f <- \(x){
     x <- matrix(x, n)
     x <- (x + t(x)) / 2
     solve(x[l1:l2, l1:l2], x[l1:l2, k1:k2]) |> t()
     }
     g <- \(x) cos(x) |> sum()
     h <- \(x) f(x) |> g()

     dput(numDeriv::grad(h, c(Sig)) + 1)
     dput(f(Sig) |> numDeriv::grad(func = g))
     */
    constexpr size_t n{9}, k1{1}, k2{2}, l1{4}, l2{7};
    constexpr double Sig[]{6.08, 1.46, -0.8, -3.43, 1.83, 4.94, -0.55, -3.13, -0.5, 1.46, 7.55, 0.54, -3.33, 0.22, 3.7, -1.85, 0.97, -0.73, -0.8, 0.54, 8.6, -0.75, -2.57, -5.09, 1.35, -0.97, 1.62, -3.43, -3.33, -0.75, 11.88, -3.46, -4.85, -0.18, -1.71, 2.55, 1.83, 0.22, -2.57, -3.46, 5.09, 3.21, 0.41, -1.04, 2.09, 4.94, 3.7, -5.09, -4.85, 3.21, 13.31, -1.73, -0.35, 0.4, -0.55, -1.85, 1.35, -0.18, 0.41, -1.73, 1.81, 0.52, 3.96, -3.13, 0.97, -0.97, -1.71, -1.04, -0.35, 0.52, 4.97, 0.6, -0.5, -0.73, 1.62, 2.55, 2.09, 0.4, 3.96, 0.6, 16.32},
                      dZ[]{-0.114426572533919, 0.536354473107216, -0.125728344982415, 0.145313181448505, 0.854131447234853, -0.753529982180174, -0.328931319719274, 0.401536257206844},
                   truth[]{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.906954228507178, 1.05943267114733, 1.33788507156456, 0.916271214399411, 1, 1, 1, 1, 1, 1.14521584164052, 0.929191677893805, 0.662213287800272, 1.10113858785874, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.906954228507178, 1.14521584164052, 1, 1.18575735592204, 0.986005257320029, 0.550830100155968, 1.15804614833241, 1, 1, 1.05943267114733, 0.929191677893805, 1, 0.986005257320029, 0.964363715691559, 1.02942700065468, 0.976126307158698, 1, 1, 1.33788507156456, 0.662213287800272, 1, 0.550830100155968, 1.02942700065468, 2.26845109143418, 0.575137856243069, 1, 1, 0.916271214399411, 1.10113858785874, 1, 1.15804614833241, 0.976126307158698, 0.575137856243069, 1.13970668414913, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
                     shift{1};

    ghqCpp::simple_mem_stack<double> mem;
    double * d_Sig{mem.get(n * n)};
    std::fill(d_Sig, d_Sig + n * n, shift);

    lp_mmcif::backprop_cond_mean
      (dZ, Sig, d_Sig, k1, k2, l1, l2, n, mem);

    for(size_t i = 0; i < n * n; ++i)
      expect_true(std::abs(d_Sig[i] - truth[i]) < std::abs(truth[i]) * 1e-8);
  }

  test_that("backprop_cond_vcov works as expected"){
    /*
     n <- 9L
     k1 <- 2L
     k2 <- 3L
     l1 <- 5L
     l2 <- 8L
     set.seed(22)
     dput(Sig <- rWishart(1, n, diag(n)) |> drop() |> round(2))
     f <- \(x){
     x <- matrix(x, n)
     x <- (x + t(x)) / 2
     x[k1:k2, k1:k2] - x[k1:k2, l1:l2] %*% solve(x[l1:l2, l1:l2], x[l1:l2, k1:k2])
     }
     g <- \(x) cos(x) |> sum()
     h <- \(x) f(x) |> g()

     dput(numDeriv::grad(h, c(Sig)) - 2)
     dput(f(Sig) |> numDeriv::grad(func = g))
     */
    constexpr size_t n{9}, k1{1}, k2{2}, l1{4}, l2{7};
    constexpr double Sig[]{6.08, 1.46, -0.8, -3.43, 1.83, 4.94, -0.55, -3.13, -0.5, 1.46, 7.55, 0.54, -3.33, 0.22, 3.7, -1.85, 0.97, -0.73, -0.8, 0.54, 8.6, -0.75, -2.57, -5.09, 1.35, -0.97, 1.62, -3.43, -3.33, -0.75, 11.88, -3.46, -4.85, -0.18, -1.71, 2.55, 1.83, 0.22, -2.57, -3.46, 5.09, 3.21, 0.41, -1.04, 2.09, 4.94, 3.7, -5.09, -4.85, 3.21, 13.31, -1.73, -0.35, 0.4, -0.55, -1.85, 1.35, -0.18, 0.41, -1.73, 1.81, 0.52, 3.96, -3.13, 0.97, -0.97, -1.71, -1.04, -0.35, 0.52, 4.97, 0.6, -0.5, -0.73, 1.62, 2.55, 2.09, 0.4, 3.96, 0.6, 16.32},
                      dZ[]{0.991989215952393, 0.0421245737377473, 0.0421245737377473, 0.990558573277914},
                   truth[]{-2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -1.00801078405628, -1.95787542626121, -2, -2.08991188571459, -2.11890910688617, -1.0202733879833, -2.31508110458744, -2, -2, -1.95787542626121, -1.00944142670963, -2, -1.44406388462062, -1.86085763871797, -2.80222728278984, -1.60482644804911, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2.08991188571459, -1.44406388462062, -2, -1.67496712785522, -1.90759364134745, -2.56650316238225, -1.74015480523204, -2, -2, -2.11890910688617, -1.86085763871797, -2, -1.90759364134745, -1.96471902209964, -2.24049467347108, -1.90265229187432, -2, -2, -1.0202733879833, -2.80222728278984, -2, -2.56650316238225, -2.24049467347108, -0.312246110311083, -2.65985185607703, -2, -2, -2.31508110458744, -1.60482644804911, -2, -1.74015480523204, -1.90265229187432, -2.65985185607703, -1.73111052666525, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
                     shift{-2};

    ghqCpp::simple_mem_stack<double> mem;
    double * d_Sig{mem.get(n * n)};
    std::fill(d_Sig, d_Sig + n * n, shift);

    lp_mmcif::backprop_cond_vcov
      (dZ, Sig, d_Sig, k1, k2, l1, l2, n, mem);

    for(size_t i = 0; i < n * n; ++i)
      expect_true(std::abs(d_Sig[i] - truth[i]) < std::abs(truth[i]) * 1e-8);
  }

  test_that("backprop_cond_vcov_rev works as expected"){
    /*
     n <- 4L
     set.seed(22)
     dput(Sig <- rWishart(1, n, diag(n)) |> drop() |> round(2))
     v1 <- rnorm(n)

     f <- \(x){
     x <- matrix(x, n)
     x <- (x + t(x)) / 2
     solve(v1 %o% v1 + solve(x))
     }
     g <- \(x) cos(x) |> sum()
     h <- \(x) f(x) |> g()

     dput(numDeriv::grad(h, c(Sig)) - 2.5)
     dput(f(Sig))
     dput(f(Sig) |> numDeriv::grad(func = g))
     */
    constexpr size_t n{4};
    constexpr double Sig[]{1.88, 0.81, -0.45, -1.28, 0.81, 2.46, 0.2, -1.03, -0.45, 0.2, 2.29, -0.31, -1.28, -1.03, -0.31, 3.29},
                       M[]{1.6865098385872, 0.384647005189191, -0.936210804617278, -0.821059137253677, 0.384647005189191, 1.52493865696598, -0.868846190128029, -0.0211018266518998, -0.936210804617278, -0.868846190128029, 1.06822756051029, 0.843247299595629, -0.821059137253678, -0.0211018266518993, 0.843247299595629, 2.20143449382433},
                      dZ[]{-0.993312658343481, -0.37523195857404, 0.805317486159614, 0.731867985509628, -0.37523195857404, -0.998948721329163, 0.763584421280266, 0.0211002580060688, 0.805317486159614, 0.763584421280266, -0.876348135694605, -0.746806641823308, 0.731867985509627, 0.0211002580060693, -0.746806641823308, -0.807651371011641},
                   truth[]{-3.29728074796648, -2.66165225279425, -1.50694505099155, -1.9139965121352, -2.66165225279425, -3.34323866627478, -1.81262168721966, -2.61217177834769, -1.50694505099155, -1.81262168721966, -4.22034621479994, -3.29296463385476, -1.9139965121352, -2.61217177834769, -3.29296463385476, -3.20766114886176},
                     shift{-2.5};

    ghqCpp::simple_mem_stack<double> mem;
    double * d_Sig{mem.get(n * n)};
    std::fill(d_Sig, d_Sig + n * n, shift);

    lp_mmcif::backprop_cond_vcov_rev(dZ, Sig, M, d_Sig, n, mem);

    for(size_t i = 0; i < n * n; ++i)
      expect_true(std::abs(d_Sig[i] - truth[i]) < std::abs(truth[i]) * 1e-8);

  }
}