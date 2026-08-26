#include "pricers/base_pricer.hpp"

#include "core/constants.hpp"
#include "math/lsmc_utils.hpp"

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

double BasePricer::price(const EuropeanOption& option) {
    const double r = option.getR();
    const double dt = option.getT() / timeSteps_;

    std::vector<std::vector<double>> assetPaths(paths_, std::vector<double>(timeSteps_ + 1));
    generateGBMPaths(option, dt, assetPaths);

    double payoffSum = 0.0;

    #pragma omp parallel for reduction(+:payoffSum)
    for (int i = 0; i < paths_; ++i) {
        payoffSum += option.payoff(assetPaths[i][timeSteps_]);
    }

    return std::exp(-r * option.getT()) * payoffSum / paths_;
}

double BasePricer::priceArithmeticAsian(const BaseOption& option) {
    const double r = option.getR();
    const double dt = option.getT() / timeSteps_;

    std::vector<std::vector<double>> assetPaths(paths_, std::vector<double>(timeSteps_ + 1));
    generateGBMPaths(option, dt, assetPaths);

    double payoffSum = 0.0;

    #pragma omp parallel for reduction(+:payoffSum)
    for (int i = 0; i < paths_; ++i) {
        double average = 0.0;

        for (int j = 1; j <= timeSteps_; ++j) {
            average += assetPaths[i][j];
        }

        average /= timeSteps_;
        payoffSum += option.payoff(average);
    }

    return std::exp(-r * option.getT()) * payoffSum / paths_;
}

double BasePricer::price(const AmericanOption& option) {
    const double K = option.getK();
    const double r = option.getR();
    const double dt = option.getT() / timeSteps_;

    std::vector<std::vector<double>> assetPaths(paths_, std::vector<double>(timeSteps_ + 1));
    std::vector<std::vector<double>> cashFlows(paths_, std::vector<double>(timeSteps_ + 1, 0.0));

    generateGBMPaths(option, dt, assetPaths);

    #pragma omp parallel for
    for (int i = 0; i < paths_; ++i) {
        cashFlows[i][timeSteps_] = option.payoff(assetPaths[i][timeSteps_]);
    }

    for (int t = timeSteps_ - 1; t >= 1; --t) {
        std::vector<double> inMoneyPaths;
        std::vector<double> inMoneyValues;
        std::vector<int> inMoneyIndices;

        for (int i = 0; i < paths_; ++i) {
            const double immediate = option.payoff(assetPaths[i][t]);

            if (immediate > 0.0) {
                inMoneyPaths.push_back(assetPaths[i][t] / K);
                inMoneyValues.push_back(cashFlows[i][t + 1] * std::exp(-r * dt));
                inMoneyIndices.push_back(i);
            }
        }

        if (inMoneyPaths.size() > OptionConstants::BASIS_FUNCTIONS_COUNT) {
            const std::vector<double> beta = MathUtils::leastSquares(inMoneyPaths, inMoneyValues);

            #pragma omp parallel for
            for (int j = 0; j < static_cast<int>(inMoneyIndices.size()); ++j) {
                const int i = inMoneyIndices[j];
                const double immediate = option.payoff(assetPaths[i][t]);
                const double continuation = MathUtils::continuationValue(assetPaths[i][t], beta, K);

                if (immediate > continuation) {
                    cashFlows[i][t] = immediate;

                    for (int s = t + 1; s <= timeSteps_; ++s) {
                        cashFlows[i][s] = 0.0;
                    }
                } else {
                    cashFlows[i][t] = cashFlows[i][t + 1];
                }
            }

            #pragma omp parallel for
            for (int i = 0; i < paths_; ++i) {
                bool inMoney = false;

                for (int index : inMoneyIndices) {
                    if (index == i) {
                        inMoney = true;
                        break;
                    }
                }

                if (!inMoney) {
                    cashFlows[i][t] = cashFlows[i][t + 1];
                }
            }
        } else {
            #pragma omp parallel for
            for (int i = 0; i < paths_; ++i) {
                cashFlows[i][t] = cashFlows[i][t + 1];
            }
        }
    }

    double optionValue = 0.0;

    #pragma omp parallel for reduction(+:optionValue)
    for (int i = 0; i < paths_; ++i) {
        optionValue += cashFlows[i][1] * std::exp(-r * dt);
    }

    return optionValue / paths_;
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
