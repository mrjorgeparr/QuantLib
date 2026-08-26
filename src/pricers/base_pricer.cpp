#include "pricers/base_pricer.hpp"

#include <cmath>
#include <omp.h>
#include <random>
#include <stdexcept>

BasePricer::BasePricer(int timeSteps, int paths)
    : timeSteps_(timeSteps),
      paths_(paths),
      gen_(std::random_device{}()) {
    if (timeSteps_ <= 0) {
        throw std::invalid_argument("Pricer time steps must be positive");
    }

    if (paths_ <= 0) {
        throw std::invalid_argument("Pricer path count must be positive");
    }
}

void BasePricer::generateGBMPaths(
    const BaseOption& option,
    double dt,
    std::vector<std::vector<double>>& assetPaths
) {
    const auto seed = gen_();
    const double S0 = option.getTicker().getS0();
    const double r = option.getR();
    const double sigma = option.getTicker().getSigma();

    #pragma omp parallel
    {
        std::mt19937 localGen(seed + omp_get_thread_num());
        std::normal_distribution<double> localNorm(0.0, 1.0);

        #pragma omp for
        for (int i = 0; i < paths_; ++i) {
            assetPaths[i][0] = S0;

            for (int j = 1; j <= timeSteps_; ++j) {
                const double dW = std::sqrt(dt) * localNorm(localGen);
                assetPaths[i][j] = assetPaths[i][j - 1]
                    * std::exp((r - 0.5 * sigma * sigma) * dt + sigma * dW);
            }
        }
    }
}
