#ifndef EUROPEAN_MONTE_CARLO_HPP
#define EUROPEAN_MONTE_CARLO_HPP

#include "instruments/european_option.hpp"

#include "pricers/base_pricer.hpp"

#include <vector>

class EuropeanPricer : public BasePricer {
private:
    double dt_;

public:
    EuropeanPricer(int timeSteps, int paths);

    double price(const EuropeanOption& option);
};

#endif // LSMC_HPP
