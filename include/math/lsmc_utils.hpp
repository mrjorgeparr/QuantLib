#ifndef MATH_UTILS_HPP
#define MATH_UTILS_HPP

#include <vector>

namespace MathUtils {
    // Basis functions for regression (Laguerre polynomials)
    std::vector<double> basisFunctions(double x);
    
    // Least squares regression solver
    std::vector<double> leastSquares(const std::vector<double>& x, 
                                   const std::vector<double>& y);
    
    // Calculate continuation value using regression coefficients
    double continuationValue(double S, const std::vector<double>& beta, double K);
    
    // Gaussian elimination solver for 4x4 system
    void gaussianElimination(std::vector<std::vector<double>>& A, 
                           std::vector<double>& b, 
                           std::vector<double>& x);
}

#endif // MATH_UTILS_HPP