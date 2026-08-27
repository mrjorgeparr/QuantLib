#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

namespace OptionConstants {
    constexpr int DEFAULT_TIME_STEPS = 50;
    constexpr int DEFAULT_MC_PATHS = 100000;
    constexpr int BASIS_FUNCTIONS_COUNT = 4;
    // Finite-difference bump sizes for the Greeks.
    //
    // The spot bumps are RELATIVE (a fraction of S0), not absolute: an absolute
    // bump of 0.01 on a spot of 100 is far below the Monte Carlo noise floor.
    // Gamma gets a larger bump than delta because the second difference divides
    // by h^2, which amplifies sampling error quadratically.
    constexpr double FINITE_DIFF_BUMP_PRICE = 0.01;      // 1% of spot  (delta)
    constexpr double FINITE_DIFF_BUMP_GAMMA = 0.05;      // 5% of spot  (gamma)
    constexpr double FINITE_DIFF_BUMP_TIME = 1.0/365.0;  // 1 day
    constexpr double FINITE_DIFF_BUMP_VOL = 0.01;        // 1 vol point
    constexpr double FINITE_DIFF_BUMP_RATE = 0.001;      // 10 basis points
    constexpr double EPSILON = 1e-10;
}

#endif // CONSTANTS_HPP

