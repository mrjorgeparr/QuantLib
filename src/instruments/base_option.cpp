#include "instruments/base_option.hpp"

#include <algorithm>
#include <stdexcept>

namespace {
void validateContract(double strike, double maturity) {
    if (strike <= 0.0) {
        throw std::invalid_argument("Option strike must be positive");
    }

    if (maturity < 0.0) {
        throw std::invalid_argument("Option maturity cannot be negative");
    }
}
}

BaseOption::BaseOption(const Ticker& ticker, double strike, double rate, double maturity, bool isPut)
    : ticker_(ticker),
      strike_(strike),
      rate_(rate),
      maturity_(maturity),
      isPut_(isPut) {
    validateContract(strike_, maturity_);
}

double BaseOption::payoff(double S) const {
    if (isPut_) {
        return std::max(strike_ - S, 0.0);
    }

    return std::max(S - strike_, 0.0);
}

void BaseOption::setTicker(const Ticker& ticker) {
    ticker_ = ticker;
}

void BaseOption::setK(double strike) {
    validateContract(strike, maturity_);
    strike_ = strike;
}

void BaseOption::setR(double rate) {
    rate_ = rate;
}

void BaseOption::setT(double maturity) {
    validateContract(strike_, maturity);
    maturity_ = maturity;
}

void BaseOption::setIsPut(bool isPut) {
    isPut_ = isPut;
}


