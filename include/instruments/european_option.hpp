#ifndef EUROPEAN_OPTION_HPP
#define EUROPEAN_OPTION_HPP

#include "instruments/base_option.hpp"

// A European option is fully described by its contract terms (BaseOption).
// Discretisation -- time steps and path count -- belongs to the pricer, not here.
// This type exists so BasePricer can overload on the exercise style.
class EuropeanOption : public BaseOption {
public:
    // Constructor
    EuropeanOption(double S0, double K, double r, double sigma, double T,
                   bool isPut = true);

    // Getters
    double getS0() const { return getTicker().getS0(); }
    double getSigma() const { return getTicker().getSigma(); }
};

#endif // EUROPEAN_OPTION_HPP
