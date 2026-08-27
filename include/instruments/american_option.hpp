#ifndef AMERICAN_OPTION_HPP
#define AMERICAN_OPTION_HPP

#include "instruments/base_option.hpp"

#include <random>

// An American option is fully described by its contract terms (BaseOption).
// Discretisation -- time steps and path count -- belongs to the pricer, not here.
// This type exists so BasePricer can overload on the exercise style.
class AmericanOption : public BaseOption {
private:
    // TODO(unused): leftover from when the instrument priced itself.
    // Path generation now lives in BasePricer::generateGBMPaths.
    std::mt19937 gen;
    std::normal_distribution<double> norm_dist;

public:
    // Constructor
    AmericanOption(double S0, double K, double r, double sigma, double T,
                   bool isPut = true);

    // Getters
    double getS0() const { return getTicker().getS0(); }
    double getSigma() const { return getTicker().getSigma(); }
};

#endif // AMERICAN_OPTION_HPP
