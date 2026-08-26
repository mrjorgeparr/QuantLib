#include "core/ticker.hpp"

Ticker::Ticker(double S0, double sigma)
    : S0(S0), sigma(sigma) {}

void Ticker::setS0(double S) { S0 = S; }
void Ticker::setSigma(double sigma) { this->sigma = sigma; }
