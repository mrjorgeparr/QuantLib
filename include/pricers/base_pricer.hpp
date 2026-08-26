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

public:
    BasePricer(int timeSteps, int paths);
    virtual ~BasePricer() = default;

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
