#include "../include/hermite.hpp"
#include "../src/overlap_binomial_expansion.cpp"
#include <iostream>
#include <cassert>
#include <cmath>
// Unit test for overlap of two 1D Hermite-Gaussians with two exponents (a, b)
// and center displacement (Qx).This test checks if build_hermite_gaussian returns that value with high accuracy.

// Utility to compare doubles
void assert_close(double val1, double val2, double tol = 1e-12) {
    assert(std::abs(val1 - val2) < tol && "Values not within tolerance!");
}

int main() {
    const double a = 0.5;  // Orbital exponent of Gaussian 'a'
    const double b = 0.5;  // Orbital exponent of Gaussian 'b'

    // === Test 1: Overlap of two s-type Gaussians (i = j = 0) ===
    {
        double Qx = 1.0;//distance between two gaussian centers
        double mu = a * b / (a + b);
        double expected = std::exp(-mu * Qx * Qx);
        double result = build_hermite_gaussian(0, 0, 0, Qx, a, b);
        std::cout << "[s] Qx = 1.0 | Expected: " << expected << ", Got: " << result << std::endl;
        assert_close(result, expected);
    }

    // === Test 2: Overlap between p-type and s-type Gaussians (i = 1, j = 0) ===
    {
        double Qx = 1.0;
        double P = a + b;
        double mu = a * b / P;
        double expected = -2.0 * mu * Qx * std::exp(-mu * Qx * Qx);  // From H1(x) = 2x
        double result = build_hermite_gaussian(1, 0, 0, Qx, a, b);
        std::cout << "[p] Qx = 1.0 | Expected: " << expected << ", Got: " << result << std::endl;
        assert_close(result, expected);
    }

    // === Test 3: Overlap between d-type and s-type Gaussians (i = 2, j = 0)  ===
    {//This one is tricky, cannot use simple closed form like in the ss and sp case
        double Qx = 1.0;
        double P = a + b;
        double expected = gaussian_overlap_1D(2, 0, a, 0, b, Qx)/std::sqrt(M_PI/P);
        double result = build_hermite_gaussian(2, 0, 0, Qx, a, b);
        std::cout << "[d] Qx = 1.0 | Expected: " << expected << ", Got: " << result << std::endl;
        assert_close(result, expected);
    }

    // === Test 4: Symmetry test for even angular momentum (i = 2, j = 0) ===
    {
        double Qx_pos = 1.0;
        double Qx_neg = -1.0;
        double result_pos = build_hermite_gaussian(2, 0, 0, Qx_pos, a, b);
        double result_neg = build_hermite_gaussian(2, 0, 0, Qx_neg, a, b);
        std::cout << "[d] Qx = +/-1.0 | +Qx: " << result_pos << ", -Qx: " << result_neg << std::endl;
        assert_close(result_pos, result_neg);
    }

    std::cout << "All tests passed." << std::endl;
    return 0;
}