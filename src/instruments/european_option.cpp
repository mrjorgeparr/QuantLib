#include "instruments/european_option.hpp"
#include "math/lsmc_utils.hpp"
#include "core/constants.hpp"
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <stdexcept>
EuropeanOption::EuropeanOption(double S0, double K, double r, double sigma, 
                                       double T, int N, int M, bool isPut)
    : BaseOption(Ticker(S0, sigma), K, r, T, isPut),
      N(N),
      M(M),
      dt(T / N) {
    if (N <= 0) {
        throw std::invalid_argument("American option time steps must be positive");
    }

    if (M <= 0) {
        throw std::invalid_argument("American option path count must be positive");
    }
}
  