#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

namespace OptionConstants {
    constexpr int DEFAULT_TIME_STEPS = 50;
    constexpr int DEFAULT_MC_PATHS = 100000;
    constexpr int BASIS_FUNCTIONS_COUNT = 4;
    constexpr double FINITE_DIFF_BUMP_PRICE = 0.01;
    constexpr double FINITE_DIFF_BUMP_TIME = 1.0/365.0;  // 1 day
    constexpr double FINITE_DIFF_BUMP_VOL = 0.01;
    constexpr double FINITE_DIFF_BUMP_RATE = 0.001;      // 10 basis points
    constexpr double EPSILON = 1e-10;
}

#endif // CONSTANTS_HPP

