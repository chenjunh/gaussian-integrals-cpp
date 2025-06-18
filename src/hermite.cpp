#include <cmath>
#include "../include/hermite.hpp"

// Recursive function to build Hermite Gaussian coefficients
//ref: https://joshuagoings.com/assets/integrals.pdf
double build_hermite_gaussian(int i,int j,int t, double Qx, double a, double b) 
{
/*
Recursive definition of Hermite Gaussian coefficients.
Returns a float.
a: orbital exponent on Gaussian 'a' (e.g. alpha in the text)
b: orbital exponent on Gaussian 'b' (e.g. beta in the text)
i,j: orbital angular momentum number on Gaussian 'a' and 'b'
t: number nodes in Hermite (depends on type of integral,
e.g. always zero for overlap integrals)
Qx: distance between origins of Gaussian 'a' and 'b'
*/ 

    double p=a+b;
    double q=a*b/p;
    if (t<0 || t>i+j) {
        return 0.0;
    }
    if (i==0 && j==0 && t==0) {
        return std::exp(-q*Qx*Qx);
    }
    
    if (j==0) {
        return 1.0/2.0/p*build_hermite_gaussian(i-1,j,t-1,Qx,a,b)
        -q*Qx/a*build_hermite_gaussian(i-1,j,t,Qx,a,b)
        +(t+1)*build_hermite_gaussian(i-1,j,t+1,Qx,a,b); 
    } 
    else {
        return 1.0/2.0/p*build_hermite_gaussian(i,j-1,t-1,Qx,a,b)
        +q*Qx/b*build_hermite_gaussian(i,j-1,t,Qx,a,b)
        +(t+1)*build_hermite_gaussian(i,j-1,t+1,Qx,a,b); 
    }
};