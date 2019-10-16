#ifndef SECANT_METHOD_H
#define SECANT_METHOD_H

#include <exception>
#include <iostream>
#include <sstream>
#include <string>



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
double secantMethod(const F& f, double x0, double x1,
        double TOL = 1.0e-12, int MAXITER = 10000)
{
    int i = 0;
    double xn = x1; 
    while (std::abs(f(xn)) > TOL && i < MAXITER)
    {
        xn = x1 - f(x1)*(x1-x0)/(f(x1)-f(x0));
        x0 = x1;
        x1 = xn;
        ++i;
    }
 
    if (i == MAXITER)
        throwSecantMethodException("Did not find root!",TOL,MAXITER);
    return xn;
}


#endif
