#include "greeks/greeks.hpp"

#include "core/constants.hpp"

#include <stdexcept>

namespace {

// Saves the common-random-number state of a pricer and restores it on scope
// exit, so computing Greeks never leaves the caller pricer reconfigured.
class SeedGuard {
public:
    explicit SeedGuard(BasePricer& pricer)
        : pricer_(pricer),
          hadSeed_(pricer.hasFixedSeed()),
          seed_(pricer.getFixedSeed()) {}

    ~SeedGuard() {
        if (hadSeed_) {
            pricer_.setSeed(seed_);
        } else {
            pricer_.clearSeed();
        }
    }

private:
    SeedGuard(const SeedGuard&);
    SeedGuard& operator=(const SeedGuard&);

    BasePricer& pricer_;
    bool hadSeed_;
    std::mt19937::result_type seed_;
};

}  // namespace

GreekBumps::GreekBumps()
    : spotRel(OptionConstants::FINITE_DIFF_BUMP_PRICE),
      gammaRel(OptionConstants::FINITE_DIFF_BUMP_GAMMA),
      vol(OptionConstants::FINITE_DIFF_BUMP_VOL),
      rate(OptionConstants::FINITE_DIFF_BUMP_RATE),
      time(OptionConstants::FINITE_DIFF_BUMP_TIME) {}

// --- scenario construction ------------------------------------------------
//
// Each of these is a copy of the contract with one input moved, built with
// BaseOption and Ticker. setT and setK run their own validation, so an
// invalid bump throws from the instrument rather than from here.

BaseOption Greeks::withSpot(double spot) const {
    BaseOption terms = option_;
    terms.setTicker(Ticker(spot, option_.getTicker().getSigma()));
    return terms;
}

BaseOption Greeks::withVol(double vol) const {
    BaseOption terms = option_;
    terms.setTicker(Ticker(option_.getTicker().getS0(), vol));
    return terms;
}

BaseOption Greeks::withRate(double rate) const {
    BaseOption terms = option_;
    terms.setR(rate);
    return terms;
}

BaseOption Greeks::withMaturity(double maturity) const {
    BaseOption terms = option_;
    terms.setT(maturity);
    return terms;
}

// Prices one leg. Pinning the seed here is what makes the two legs of every
// difference share their random draws.
double Greeks::priceOf(const BaseOption& terms) {
    if (useCommonRandomNumbers_) {
        pricer_.setSeed(seed_);
    }

    return reprice_(terms);
}

double Greeks::spotBump() const {
    if (spot() <= 0.0) {
        throw std::invalid_argument("Greeks require a positive spot");
    }

    return bumps_.spotRel * spot();
}

double Greeks::gammaBump() const {
    if (spot() <= 0.0) {
        throw std::invalid_argument("Greeks require a positive spot");
    }

    return bumps_.gammaRel * spot();
}

// --- the Greeks -----------------------------------------------------------

double Greeks::price() {
    SeedGuard guard(pricer_);
    return priceOf(option_);
}

// dV/dS, central difference.
double Greeks::delta() {
    SeedGuard guard(pricer_);

    const double h = spotBump();

    const double up = priceOf(withSpot(spot() + h));
    const double down = priceOf(withSpot(spot() - h));

    return (up - down) / (2.0 * h);
}

// d2V/dS2, symmetric second difference.
double Greeks::gamma() {
    SeedGuard guard(pricer_);

    const double h = gammaBump();

    const double up = priceOf(withSpot(spot() + h));
    const double mid = priceOf(option_);
    const double down = priceOf(withSpot(spot() - h));

    return (up - 2.0 * mid + down) / (h * h);
}

// dV/dsigma, per 1.0 of volatility. One-sided when sigma sits below the bump.
double Greeks::vega() {
    SeedGuard guard(pricer_);

    const double h = bumps_.vol;
    const double up = priceOf(withVol(vol() + h));

    if (vol() <= h) {
        return (up - priceOf(option_)) / h;
    }

    const double down = priceOf(withVol(vol() - h));

    return (up - down) / (2.0 * h);
}

// dV/dr, per 1.0 of rate.
double Greeks::rho() {
    SeedGuard guard(pricer_);

    const double h = bumps_.rate;
    const double r = option_.getR();

    const double up = priceOf(withRate(r + h));
    const double down = priceOf(withRate(r - h));

    return (up - down) / (2.0 * h);
}

// Theta is the derivative with respect to CALENDAR time, and time to maturity
// runs the other way: theta = dV/dt = -dV/dT. Negative for most long options.
// Near expiry, T - h would go negative -- and BaseOption::setT would throw --
// so fall back to a one-sided difference.
double Greeks::theta() {
    SeedGuard guard(pricer_);

    const double h = bumps_.time;
    const double T = option_.getT();
    const double up = priceOf(withMaturity(T + h));

    if (T <= h) {
        return -(up - priceOf(option_)) / h;
    }

    const double down = priceOf(withMaturity(T - h));

    return -(up - down) / (2.0 * h);
}

// Computes everything in one sweep, reusing the base price for gamma and for
// the one-sided fallbacks: ten pricer calls instead of the twelve the
// individual accessors would take.
GreekSet Greeks::all() {
    SeedGuard guard(pricer_);

    GreekSet greeks;
    greeks.price = priceOf(option_);

    const double S = spot();
    const double sigma = vol();
    const double r = option_.getR();
    const double T = option_.getT();

    const double hSpot = spotBump();
    greeks.delta = (priceOf(withSpot(S + hSpot))
                    - priceOf(withSpot(S - hSpot))) / (2.0 * hSpot);

    const double hGamma = gammaBump();
    greeks.gamma = (priceOf(withSpot(S + hGamma))
                    - 2.0 * greeks.price
                    + priceOf(withSpot(S - hGamma))) / (hGamma * hGamma);

    const double hVol = bumps_.vol;
    const double volUp = priceOf(withVol(sigma + hVol));
    if (sigma <= hVol) {
        greeks.vega = (volUp - greeks.price) / hVol;
    } else {
        greeks.vega = (volUp - priceOf(withVol(sigma - hVol))) / (2.0 * hVol);
    }

    const double hRate = bumps_.rate;
    greeks.rho = (priceOf(withRate(r + hRate))
                  - priceOf(withRate(r - hRate))) / (2.0 * hRate);

    const double hTime = bumps_.time;
    const double timeUp = priceOf(withMaturity(T + hTime));
    if (T <= hTime) {
        greeks.theta = -(timeUp - greeks.price) / hTime;
    } else {
        greeks.theta = -(timeUp - priceOf(withMaturity(T - hTime))) / (2.0 * hTime);
    }

    return greeks;
}
