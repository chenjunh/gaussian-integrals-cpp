#include "../src/overlap_binomial_expansion.cpp"
#include "../src/overlap_hermite_recursive.cpp"
#include <iostream>
#include <cassert>
#include <cmath>
// Utility to compare doubles
void assert_close(double val1, double val2, double tol = 1e-12) {
    assert(std::abs(val1 - val2) < tol && "Values not within tolerance!");
}

int main() {
    const double a = 0.5;  // Orbital exponent of Gaussian 'a'
    const double b = 0.5;  // Orbital exponent of Gaussian 'b'
    const std::vector<double> A = {0,0,0};//Orbital center A
    const std::vector<double> B = {1.0,1.0,1.0};//Orbital center B
    
    // === Test 1: compare Overlap of two s-type Gaussians (i = j = 0) ===
    {
        std::vector<int> lmn1 = {0,0,0};//Angular momentum of orbital 1
        std::vector<int> lmn2 = {0,0,0};//Angular momentum of orbital 2
        double expected = gaussian_overlap_1D(lmn1[0], lmn2[0], a, A[0], b, B[0]) 
                        *gaussian_overlap_1D(lmn1[1], lmn2[1], a, A[1], b, B[1])
                        *gaussian_overlap_1D(lmn1[2], lmn2[2], a, A[2], b, B[2]);
        double result = goverlap(a, lmn1, A, b, lmn2, B);
        std::cout << "[ss] | Expected: " << expected << ", Got: " << result << std::endl;
        assert_close(result, expected);
    }

    // === Test 2: Overlap between p-type and s-type Gaussians (i = 1, j = 0) ===
    {
        std::vector<int> lmn1 = {0,0,0};//s
        std::vector<int> lmn2 = {1,0,0};//px
        double expected = gaussian_overlap_1D(lmn1[0], lmn2[0], a, A[0], b, B[0]) 
                        *gaussian_overlap_1D(lmn1[1], lmn2[1], a, A[1], b, B[1])
                        *gaussian_overlap_1D(lmn1[2], lmn2[2], a, A[2], b, B[2]);
        double result = goverlap(a, lmn1, A, b, lmn2, B);
        std::cout << "[sp] | Expected: " << expected << ", Got: " << result << std::endl;
        assert_close(result, expected);
    }

        // === Test 3: Overlap between d-type and s-type Gaussians (i = 2, j = 0) ===
    {
        std::vector<int> lmn1 = {0,0,0};//s
        std::vector<int> lmn2 = {2,0,0};//dx^2
        double expected = gaussian_overlap_1D(lmn1[0], lmn2[0], a, A[0], b, B[0]) 
                        *gaussian_overlap_1D(lmn1[1], lmn2[1], a, A[1], b, B[1])
                        *gaussian_overlap_1D(lmn1[2], lmn2[2], a, A[2], b, B[2]);
        double result = goverlap(a, lmn1, A, b, lmn2, B);
        std::cout << "[sd] | Expected: " << expected << ", Got: " << result << std::endl;
        assert_close(result, expected);
    }

    std::cout << "All tests passed." << std::endl;
    return 0;
}