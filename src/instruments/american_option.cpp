#include "instruments/american_option.hpp"

AmericanOption::AmericanOption(double S0, double K, double r, double sigma,
                               double T, bool isPut)
    : BaseOption(Ticker(S0, sigma), K, r, T, isPut),
      gen(std::random_device{}()),
      norm_dist(0.0, 1.0) {}
