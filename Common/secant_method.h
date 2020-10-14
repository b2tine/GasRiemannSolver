#ifndef SECANT_METHOD_H
#define SECANT_METHOD_H

#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>


class NoConvergenceException : public std::runtime_error
{
    public:

        NoConvergenceException(const std::string& arg, double TOL, int MAXITER)
            : std::runtime_error{arg}
        {
            std::ostringstream oss;
            oss << "ERROR: method did not converge after " << MAXITER
                << " iterations using tolerance " << TOL << "\n";
            message = oss.str();
        }

        NoConvergenceException(const std::string& arg, double TOL, int MAXITER,
                const std::string& function, const std::string& file, int line)
            : NoConvergenceException{arg,TOL,MAXITER}
        {
            std::ostringstream oss;
            oss << function << "() : " << file
                << " : " << "line " << line << " : " << arg << "\n";
            message = oss.str() + message;
        }

        const char* what() const noexcept override
        {
            return message.c_str();
        }
    
    private:
    
        std::string message;
    
};

#define throwSecantMethodException(arg,tol,maxiter) \
    throw NoConvergenceException(arg,tol,maxiter,__FUNCTION__,__FILE__,__LINE__)

template<typename F>
double secantMethod(const F& f, double xa, double xb,
        double TOL = 1.0e-07, int MAXITER = 20000)
{
    /*
    double x0,x1;
    if (fabs(xb - xa) <= TOL)
    {
        x0 = 0.001;
        x1 = 100.0;
        
    }*/
    
    double x0 = xa;
    double x1 = xb;
    if (fabs(xb - xa) <= TOL)
    {
        x0 = 0.5*(xb + xa);
        x1 = 0.5*(xb + xa) + 0.25*(xa - xb + 1.0e-06);
    }

    double xn = x0;
    for (int i = 0; i < MAXITER; ++i)
    {
        xn = x1 - f(x1)*(x1-x0)/(f(x1)-f(x0));
        x0 = x1;
        x1 = xn;
        ++i;

        double feval = f(xn);
        if (fabs(feval) < TOL)
            return xn;;
    }
    throwSecantMethodException("Did not find root!",TOL,MAXITER);
}


#endif
