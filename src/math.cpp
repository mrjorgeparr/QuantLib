#include "math/lsmc_utils.hpp"
#include "core/constants.hpp"
#include <cmath>
#include <algorithm>

namespace MathUtils {
    
    std::vector<double> basisFunctions(double x) {
        std::vector<double> basis(OptionConstants::BASIS_FUNCTIONS_COUNT);
        basis[0] = 1.0;
        basis[1] = std::exp(-x/2.0);
        basis[2] = std::exp(-x/2.0) * (1.0 - x);
        basis[3] = std::exp(-x/2.0) * (1.0 - 2.0*x + x*x/2.0);
        return basis;
    }
    
    double continuationValue(double S, const std::vector<double>& beta, double K) {
        std::vector<double> basis = basisFunctions(S/K);
        double value = 0.0;
        for (int i = 0; i < OptionConstants::BASIS_FUNCTIONS_COUNT; i++) {
            value += beta[i] * basis[i];
        }
        return value;
    }
    
    void gaussianElimination(std::vector<std::vector<double>>& A, 
                           std::vector<double>& b, 
                           std::vector<double>& x) {
        int n = static_cast<int>(b.size());
        
        // Forward elimination
        for (int i = 0; i < n; i++) {
            // Find pivot
            int maxRow = i;
            for (int k = i + 1; k < n; k++) {
                if (std::abs(A[k][i]) > std::abs(A[maxRow][i])) {
                    maxRow = k;
                }
            }
            std::swap(A[i], A[maxRow]);
            std::swap(b[i], b[maxRow]);
            
            // Make all rows below this one 0 in current column
            for (int k = i + 1; k < n; k++) {
                if (std::abs(A[i][i]) < OptionConstants::EPSILON) continue;
                double c = A[k][i] / A[i][i];
                for (int j = i; j < n; j++) {
                    if (i == j) {
                        A[k][j] = 0;
                    } else {
                        A[k][j] -= c * A[i][j];
                    }
                }
                b[k] -= c * b[i];
            }
        }
        
        // Back substitution
        for (int i = n - 1; i >= 0; i--) {
            x[i] = b[i];
            for (int j = i + 1; j < n; j++) {
                x[i] -= A[i][j] * x[j];
            }
            if (std::abs(A[i][i]) > OptionConstants::EPSILON) {
                x[i] /= A[i][i];
            }
        }
    }
    
    std::vector<double> leastSquares(const std::vector<double>& x, 
                                   const std::vector<double>& y) {
        int n = static_cast<int>(x.size());
        if (n == 0) return std::vector<double>(OptionConstants::BASIS_FUNCTIONS_COUNT, 0.0);
        
        const int p = OptionConstants::BASIS_FUNCTIONS_COUNT;
        
        // Build design matrix X
        std::vector<std::vector<double>> X(n, std::vector<double>(p));
        for (int i = 0; i < n; i++) {
            std::vector<double> basis = basisFunctions(x[i]);
            for (int j = 0; j < p; j++) {
                X[i][j] = basis[j];
            }
        }
        
        // Build X^T * X
        std::vector<std::vector<double>> XTX(p, std::vector<double>(p, 0.0));
        for (int i = 0; i < p; i++) {
            for (int j = 0; j < p; j++) {
                for (int k = 0; k < n; k++) {
                    XTX[i][j] += X[k][i] * X[k][j];
                }
            }
        }
        
        // Build X^T * y
        std::vector<double> XTy(p, 0.0);
        for (int i = 0; i < p; i++) {
            for (int k = 0; k < n; k++) {
                XTy[i] += X[k][i] * y[k];
            }
        }
        
        // Solve using Gaussian elimination
        std::vector<double> beta(p, 0.0);
        gaussianElimination(XTX, XTy, beta);
        
        return beta;
    }
}
