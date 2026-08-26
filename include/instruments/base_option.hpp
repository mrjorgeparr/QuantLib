#ifndef BASE_OPTION_HPP
#define BASE_OPTION_HPP

#include "core/ticker.hpp"



class BaseOption {
private:
    Ticker ticker_;
    double strike_;
    double rate_;
    double maturity_;
    bool isPut_;
    

public:
    BaseOption(const Ticker& ticker, double strike, double rate, double maturity, bool isPut);
    virtual ~BaseOption() = default;
    
    double getK() const { return strike_; }
    double getR() const { return rate_; }
    double getT() const { return maturity_; }
    const Ticker& getTicker() const { return ticker_; }
    bool getIsPut() const { return isPut_; }

    void setTicker(const Ticker& ticker);
    void setK(double strike);
    void setR(double rate);
    void setT(double maturity);
    void setIsPut(bool isPut);

    double payoff(double S) const;

    
};

#endif // BASE_OPTION_HPP
