#ifndef BASE_PRICER_HPP
#define BASE_PRICER_HPP

#include "instruments/american_option.hpp"
#include "instruments/base_option.hpp"
#include "instruments/european_option.hpp"

#include <random>
#include <vector>

class BasePricer {
protected:
    int timeSteps_;
    int paths_;

    std::mt19937 gen_;

    // When set, every pricing reuses identical random draws. Finite-difference
    // Greeks require this (common random numbers): without it the sampling
    // error of the two legs does not cancel and swamps the derivative.
    std::mt19937::result_type fixedSeed_;
    bool useFixedSeed_;

public:
    BasePricer(int timeSteps, int paths);
    virtual ~BasePricer() = default;

    // Common-random-number control.
    void setSeed(std::mt19937::result_type seed);
    void clearSeed();
    bool hasFixedSeed() const { return useFixedSeed_; }
    std::mt19937::result_type getFixedSeed() const { return fixedSeed_; }

    double price(const EuropeanOption& option);
    double price(const AmericanOption& option);
    double priceArithmeticAsian(const BaseOption& option);

protected:
    void generateGBMPaths(
        const BaseOption& option,
        double dt,
        std::vector<std::vector<double>>& assetPaths
    );
};

#endif // BASE_PRICER_HPP
