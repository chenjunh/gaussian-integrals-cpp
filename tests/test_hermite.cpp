#include "../include/hermite.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

int main() {
    double a = 0.5, b = 0.5, Qx = 0.0;
    double expected = std::exp(-a * b / (a + b) * Qx * Qx);
    double result = build_hermite_gaussian(0, 0, 0, Qx, a, b);

    std::cout << "Expected: " << expected << ", Got: " << result << std::endl;
    assert(std::abs(result - expected) < 1e-12);
    return 0;
}
