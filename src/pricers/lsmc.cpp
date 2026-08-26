#include "pricers/lsmc.hpp"

#include "core/constants.hpp"
#include "math/lsmc_utils.hpp"

#include <cmath>
#include <omp.h>

LSMC::LSMC(int timeSteps, int paths)
    : BasePricer(timeSteps, paths),
      dt_(0.0) {}

double LSMC::price(const AmericanOption& option) {
    const double K = option.getK();
    const double r = option.getR();
    dt_ = option.getT() / timeSteps_;

    std::vector<std::vector<double>> assetPaths(paths_, std::vector<double>(timeSteps_ + 1));
    std::vector<std::vector<double>> cashFlows(paths_, std::vector<double>(timeSteps_ + 1, 0.0));

    generateGBMPaths(option, dt_, assetPaths);

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
                inMoneyValues.push_back(cashFlows[i][t + 1] * std::exp(-r * dt_));
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
        optionValue += cashFlows[i][1] * std::exp(-r * dt_);
    }

    return optionValue / paths_;
}
