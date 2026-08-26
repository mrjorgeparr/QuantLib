#ifndef AMERICAN_OPTION_HPP
#define AMERICAN_OPTION_HPP

#include "instruments/base_option.hpp"

#include <vector>
#include <random>

class AmericanOption : public BaseOption {
private:
    int N;              // Number of time steps
    int M;              // Number of Monte Carlo paths
    double dt;          // Time step
    
    // Random number generation
    std::mt19937 gen;
    std::normal_distribution<double> norm_dist;
    
public:
    // Constructor
    AmericanOption(double S0, double K, double r, double sigma, double T, 
                       int N, int M, bool isPut = true);
    

    // Getters
    double getS0() const { return getTicker().getS0(); }
    double getSigma() const { return getTicker().getSigma(); }
    

};

#endif // AMERICAN_OPTION_HPP
