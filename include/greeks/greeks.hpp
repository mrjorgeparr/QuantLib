#ifndef GREEKS_HPP
#define GREEKS_HPP

#include "core/constants.hpp"
#include "core/ticker.hpp"
#include "instruments/base_option.hpp"
#include "pricers/base_pricer.hpp"

#include <functional>
#include <random>

// Every sensitivity Greeks computes, in one aggregate.
//
// Units: delta and gamma are per unit of spot. vega is per 1.0 of volatility
// (divide by 100 for the usual "per vol point"). rho is per 1.0 of rate
// (divide by 10000 for "per basis point"). theta is per year (divide by 365
// for "per day").
struct GreekSet {
    double price;
    double delta;
    double gamma;
    double vega;
    double rho;
    double theta;
};

// Finite-difference step sizes. The spot bumps are fractions of S0; the rest
// are absolute. Defaults come from OptionConstants.
struct GreekBumps {
    double spotRel;   // delta
    double gammaRel;  // gamma (larger: the second difference divides by h^2)
    double vol;
    double rate;
    double time;

    GreekBumps();
};

// Finite-difference Greeks for any instrument a pricer can price.
//
// The contract lives in a BaseOption member, and bumped scenarios are built
// with BaseOption and Ticker themselves -- setTicker, setR, setT, and their
// validation. Nothing about the contract is re-modelled here.
//
// The one thing BaseOption cannot do is choose a pricer overload.
// BasePricer::price is resolved at compile time from the STATIC type of its
// argument, so a BaseOption alone matches nothing. The templated constructor
// therefore captures the concrete instrument once, in a std::function that
// transplants bumped terms onto a correctly-typed copy and calls price(). That
// function is the entire cost of the overload set; everything else in this
// class is ordinary code in greeks.cpp.
//
// Accuracy: every Greek is a difference of two Monte Carlo estimates divided
// by a small bump. Unless both legs use the SAME random draws the sampling
// error does not cancel and the division amplifies what is left -- for delta
// the noise is ~1400x larger without common random numbers, enough to get the
// sign wrong. CRN is therefore ON by default.
//
// Lifetime: holds a reference to the pricer, and the captured lambda refers to
// it too. The pricer must outlive this object.
class Greeks {
public:
    template <typename OptionT>
    Greeks(const OptionT& option,
           BasePricer& pricer,
           const GreekBumps& bumps = GreekBumps());

    // --- the Greeks -------------------------------------------------------
    // Cost in pricer calls: price 1, delta 2, gamma 3, vega 2, rho 2,
    // theta 2. all() reuses the base price throughout: 10 calls, not 12.

    double price();
    double delta();
    double gamma();
    double vega();
    double rho();
    double theta();
    GreekSet all();

    // --- configuration ----------------------------------------------------

    // The contract being differentiated. Reading it goes straight through
    // BaseOption; there is no second copy of the terms in this class.
    const BaseOption& option() const { return option_; }

    // Replaces the contract terms. The instrument type captured at
    // construction is unchanged -- only the terms move.
    void setOption(const BaseOption& option) { option_ = option; }

    const GreekBumps& bumps() const { return bumps_; }
    void setBumps(const GreekBumps& bumps) { bumps_ = bumps; }

    bool commonRandomNumbers() const { return useCommonRandomNumbers_; }
    void setCommonRandomNumbers(bool on) { useCommonRandomNumbers_ = on; }

    std::mt19937::result_type seed() const { return seed_; }
    void setSeed(std::mt19937::result_type seed) { seed_ = seed; }

protected:
    BaseOption option_;
    BasePricer& pricer_;
    GreekBumps bumps_;
    bool useCommonRandomNumbers_;
    std::mt19937::result_type seed_;

    // Prices a bumped set of terms on the instrument captured at construction.
    // This is the only place the pricer overload set is reached.
    std::function<double(const BaseOption&)> reprice_;

    // Scenario construction, built out of BaseOption and Ticker.
    BaseOption withSpot(double spot) const;
    BaseOption withVol(double vol) const;
    BaseOption withRate(double rate) const;
    BaseOption withMaturity(double maturity) const;

    double priceOf(const BaseOption& terms);

    double spot() const { return option_.getTicker().getS0(); }
    double vol() const { return option_.getTicker().getSigma(); }

    double spotBump() const;
    double gammaBump() const;

private:
    template <typename OptionT>
    void capture(const OptionT& option, BasePricer& pricer);
};

// --- constructors ---------------------------------------------------------
//
// The only templated code in this class: enough to remember which price()
// overload to call, and nothing more.

template <typename OptionT>
void Greeks::capture(const OptionT& option, BasePricer& pricer) {
    BasePricer* target = &pricer;

    // `option` is captured at its exact type, so target->price(concrete) below
    // resolves to the matching overload. Copying it per call also preserves
    // any instrument-specific state a future contract might carry, which a
    // sliced BaseOption would lose.
    reprice_ = [option, target](const BaseOption& terms) {
        OptionT concrete = option;
        concrete.setTicker(terms.getTicker());
        concrete.setK(terms.getK());
        concrete.setR(terms.getR());
        concrete.setT(terms.getT());
        concrete.setIsPut(terms.getIsPut());
        return target->price(concrete);
    };
}

template <typename OptionT>
Greeks::Greeks(const OptionT& option, BasePricer& pricer, const GreekBumps& bumps)
    : option_(option),
      pricer_(pricer),
      bumps_(bumps),
      useCommonRandomNumbers_(true),
      seed_(std::random_device{}()) {
    capture(option, pricer);
}

#endif // GREEKS_HPP
