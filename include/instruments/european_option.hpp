#ifndef EUROPEAN_OPTION_HPP
#define EUROPEAN_OPTION_HPP

#include "instruments/base_option.hpp"

//#include <vector>
//#include <random>

class EuropeanOption : public BaseOption {
private:
    int N;              // Number of time steps
    int M;              // Number of Monte Carlo paths
    double dt;          // Time step
    

    
public:
    // Constructor
    EuropeanOption(double S0, double K, double r, double sigma, double T, 
                       int N, int M, bool isPut = true);
    

    // Getters
    double getS0() const { return getTicker().getS0(); }
    double getSigma() const { return getTicker().getSigma(); }
    

};

#endif // AMERICAN_OPTION_HPP
