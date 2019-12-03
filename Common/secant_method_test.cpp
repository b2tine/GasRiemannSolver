#include "secant_method.h"


double polynomialFunction(double x)
{
    return x*x - 1.0;
}

struct CubicPolynomial
{
    double A, B, C, D;

    CubicPolynomial(double a, double b, double c, double d)
        : A{a}, B{b}, C{c}, D{d}
    {}

    double operator () (double x) const
    {
        return A*x*x*x + B*x*x + C*x + D;
    }
};



int main()
{
    CubicPolynomial f(0,1,0,-1);
    std::cout << "f(1) = " << f(1) << "\n";

    double root = secantMethod(f,0,0.5); //equivalent to instantiation on next line
    double rootv = secantMethod<CubicPolynomial>(f,0,0.5);
    double proot = secantMethod(polynomialFunction,0,0.5);
    
    std::cout << "root = " << root << "\n";
    std::cout << "rootv = " << rootv << "\n";
    std::cout << "proot = " << proot << "\n";

    return 0;
}



