#include <iostream>
#include <vector>
#include <cmath>
#include <limits>

// Function to compute partial derivatives
void computePartialDerivatives(const std::vector<std::vector<double>>& implied_vol_surface,
                               const std::vector<double>& strikes,
                               const std::vector<double>& maturities,
                               std::vector<std::vector<double>>& d_sigma_dT,
                               std::vector<std::vector<double>>& d_sigma_dK,
                               std::vector<std::vector<double>>& d2_sigma_dK2) {
    int m = implied_vol_surface.size();    // number of maturities
    int n = implied_vol_surface[0].size(); // number of strikes

    // Compute d_sigma_dT using central difference
    for (int i = 1; i < m - 1; ++i) {
        for (int j = 0; j < n; ++j) {
            d_sigma_dT[i][j] = (implied_vol_surface[i + 1][j] - implied_vol_surface[i - 1][j]) /
                               (maturities[i + 1] - maturities[i - 1]);
        }
    }

    // Compute d_sigma_dK using central difference
    for (int i = 0; i < m; ++i) {
        for (int j = 1; j < n - 1; ++j) {
            d_sigma_dK[i][j] = (implied_vol_surface[i][j + 1] - implied_vol_surface[i][j - 1]) /
                               (strikes[j + 1] - strikes[j - 1]);
        }
    }

    // Compute d2_sigma_dK2 using central difference
    for (int i = 0; i < m; ++i) {
        for (int j = 1; j < n - 1; ++j) {
            d2_sigma_dK2[i][j] = (implied_vol_surface[i][j + 1] - 2 * implied_vol_surface[i][j] + implied_vol_surface[i][j - 1]) /
                                 std::pow(strikes[j + 1] - strikes[j - 1], 2);
        }
    }
}

// Function to calculate local volatility
std::vector<std::vector<double>> calculateLocalVolatility(
    const std::vector<std::vector<double>>& implied_vol_surface,
    const std::vector<std::vector<double>>& d_sigma_dT,
    const std::vector<std::vector<double>>& d_sigma_dK,
    const std::vector<std::vector<double>>& d2_sigma_dK2,
    const std::vector<double>& strikes,
    const std::vector<double>& maturities,
    double r) {

    int m = implied_vol_surface.size();
    int n = implied_vol_surface[0].size();
    std::vector<std::vector<double>> local_vol_surface(m, std::vector<double>(n, 0.0));

    for (int i = 1; i < m - 1; ++i) {
        for (int j = 1; j < n - 1; ++j) {
            double K = strikes[j];
            double T = maturities[i];
            double sigma_imp = implied_vol_surface[i][j];
            double term1 = d_sigma_dT[i][j];
            double term2 = r * K * d_sigma_dK[i][j];
            double term3 = 0.5 * K * K * d2_sigma_dK2[i][j];
            if (term3 != 0) {
                local_vol_surface[i][j] = std::sqrt((term1 + term2) / term3);
            } else {
                local_vol_surface[i][j] = std::numeric_limits<double>::quiet_NaN(); // Not-a-Number for undefined values
            }
        }
    }

    return local_vol_surface;
}

int main() {
    // Example implied volatility surface (dummy data)
    std::vector<double> strikes = {90, 100, 110};
    std::vector<double> maturities = {0.1, 0.5, 1.0};
    std::vector<std::vector<double>> implied_vols = {
        {0.2, 0.25, 0.3},
        {0.18, 0.22, 0.28},
        {0.16, 0.2, 0.25}
    };

    // Create matrices for partial derivatives
    std::vector<std::vector<double>> d_sigma_dT(maturities.size(), std::vector<double>(strikes.size(), 0.0));
    std::vector<std::vector<double>> d_sigma_dK(maturities.size(), std::vector<double>(strikes.size(), 0.0));
    std::vector<std::vector<double>> d2_sigma_dK2(maturities.size(), std::vector<double>(strikes.size(), 0.0));

    // Compute partial derivatives
    computePartialDerivatives(implied_vols, strikes, maturities, d_sigma_dT, d_sigma_dK, d2_sigma_dK2);

    // Define the risk-free rate
    double r = 0.05;

    // Calculate local volatility
    std::vector<std::vector<double>> local_vol_surface = calculateLocalVolatility(
        implied_vols, d_sigma_dT, d_sigma_dK, d2_sigma_dK2, strikes, maturities, r);

    // Print local volatility surface
    std::cout << "Local Volatility Surface:" << std::endl;
    for (const auto& row : local_vol_surface) {
        for (double value : row) {
            if (std::isnan(value)) {
                std::cout << "NaN ";
            } else {
                std::cout << value << " ";
            }
        }
        std::cout << std::endl;
    }

    return 0;
}
