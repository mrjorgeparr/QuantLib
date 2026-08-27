#include "pricers/european_pricer.hpp"

#include "core/constants.hpp"
#include "math/lsmc_utils.hpp"

#include <cmath>
#include <omp.h>
#include <random>

EuropeanPricer::EuropeanPricer(int timeSteps, int paths)
    : BasePricer(timeSteps, paths),
      dt_(0.0) {}

double EuropeanPricer::price_mc(const EuropeanOption& option) {
    const auto seed = gen_();
    const double r = option.getR();
    dt_ = option.getT() / timeSteps_;
    double optionValue = 0.0;

    // DEPRECATED: STORING WHOLE PATHS IS UNNECESARY
    /*std::vector<std::vector<double>> assetPaths(paths_, std::vector<double>(timeSteps_ + 1));
    generateGBMPaths(option, dt_, assetPaths);*/

    std::vector<double> payoffs(paths_);  
    // aggregate the values of the threads into the optionValue safely
    #pragma omp parallel for 

    // 
    std::mt19937 localGen(seed + omp_get_thread_num());
    std::normal_distribution<double> localNorm(0.0, 1.0);
    for (int i = 0; i < paths_; ++i) {

        // instantiate local thread gen key
        
        
        const double S = option.getTicker().getS0();
        // initiate starting running value to S0
         
        for (int  j = 0; i < paths_; ++j){
            //reuse payoffs
            const double dW = std::sqrt(dt_) * localNorm(localGen);
            payoffs[j] = payoffs

        }
        // take value at expiry
        optionValue += assetPaths[i][paths_];
    }

    // apply discount and compute expectation
    return optionValue * std::exp(-r * dt_) / paths_;
}
