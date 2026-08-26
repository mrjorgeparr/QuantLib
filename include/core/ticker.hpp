#ifndef TICKER_HPP
#define TICKER_HPP


class Ticker {
private:
    double S0;          // Initial stock price
    double sigma;       // Volatility is asset related

public:
    // Constructor
    Ticker(double S0, double sigma);

    // Getters
    double getS0() const { return S0; }
    double getSigma() const { return sigma; }

    // setters
    void setS0(double S);
    void setSigma(double sigma);
};

#endif
