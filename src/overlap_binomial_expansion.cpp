#include <iostream>
#include <cmath>

// Compute double factorial of even integer n (n!! = n * (n-2) * ... 2 or 1)
unsigned long long double_factorial(int n) {
    if (n <= 0) return 1;
    unsigned long long result = 1;
    for (int k = n; k > 1; k -= 2)
        result *= k;
    return result;
}

// Binomial coefficient (n choose k), naive implementation
unsigned long long binomial(int n, int k) {
    if (k > n) return 0;
    if (k == 0 || k == n) return 1;
    unsigned long long res = 1;
    for (int i = 1; i <= k; ++i) {
        res = res * (n - (k - i)) / i;
    }
    return res;
}

// Moment integral I_n = ∫ (x-P)^n e^{-γ (x-P)^2} dx over (-inf, inf)
double moment_integral(int n, double gamma) {
    if (n % 2 != 0) return 0.0; // odd moment = 0
    // even moments
    int m = n / 2;
    unsigned long long df = double_factorial(n - 1);
    double val = df * std::sqrt(M_PI) / (std::pow(2.0, m) * std::pow(gamma, m + 0.5));
    return val;
}

// Overlap integral S_{l,l'} for 1D GTFs with angular momentum l and l'
double gaussian_overlap_1D(int l, int lp,
                           double alpha_a, double A,
                           double alpha_b, double B) {
    double gamma = alpha_a + alpha_b;
    double P = (alpha_a * A + alpha_b * B) / gamma;
    double diffAB = A - B;
    double pre_exp = std::exp(- (alpha_a * alpha_b / gamma) * diffAB * diffAB);

    double S = 0.0;

    for (int i = 0; i <= l; ++i) {
        unsigned long long bin_l_i = binomial(l, i);
        double PA_pow = std::pow(P - A, l - i);
        for (int j = 0; j <= lp; ++j) {
            unsigned long long bin_lp_j = binomial(lp, j);
            double PB_pow = std::pow(P - B, lp - j);

            int n = i + j; // total power of (x-P)^n
            double mom = moment_integral(n, gamma);

            S += bin_l_i * bin_lp_j * PA_pow * PB_pow * mom;
        }
    }

    return pre_exp * S;
}
