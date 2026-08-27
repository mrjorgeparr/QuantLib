#include "instruments/european_option.hpp"

EuropeanOption::EuropeanOption(double S0, double K, double r, double sigma,
                               double T, bool isPut)
    : BaseOption(Ticker(S0, sigma), K, r, T, isPut) {}
