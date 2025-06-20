#include <cmath>
#include "../include/hermite.hpp"
#include <vector>

double goverlap(double a, std::vector<int> lmn1, std::vector<double> A, double b, std::vector<int> lmn2, std::vector<double> B) {
    /*return gaussian overlap*/
    /*Evaluates overlap integral between two Gaussians
    Returns a float.
    a: orbital exponent on Gaussian 'a' (e.g. alpha in the text)
    b: orbital exponent on Gaussian 'b' (e.g. beta in the text)
    lmn1: int vector containing orbital angular momentum (e.g. (1,0,0))
    for Gaussian 'a'
    lmn2: int vector containing orbital angular momentum for Gaussian 'b'
    A: list containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]
    B: list containing origin of Gaussian 'b'*/
    unsigned long long int l1=lmn1[0];//shell angular momentum on Gaussian 'a'
    unsigned long long int m1=lmn1[1];
    unsigned long long int n1=lmn1[2];
    unsigned long long int l2=lmn2[0];//shell angular momentum on Gaussian 'b'
    unsigned long long int m2=lmn2[1];
    unsigned long long int n2=lmn2[2];
    double S1=build_hermite_gaussian(l1,l2,0,A[0]-B[0],a,b);//X
    double S2=build_hermite_gaussian(m1,m2,0,A[1]-B[1],a,b);//Y
    double S3=build_hermite_gaussian(n1,n2,0,A[2]-B[2],a,b);//Z
    double S=S1*S2*S3*std::pow(std::sqrt(M_PI)/(a+b),3);
    return S;
};