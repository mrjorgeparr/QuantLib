#ifndef LSMC_HPP
#define LSMC_HPP

#include "instruments/american_option.hpp"
#include "pricers/base_pricer.hpp"

#include <vector>

class LSMC : public BasePricer {
private:
    double dt_;

public:
    LSMC(int timeSteps, int paths);

    double price(const AmericanOption& option);
};

#endif // LSMC_HPP
