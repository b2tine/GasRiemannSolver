#ifndef RIEMANN_PROBLEM_H
#define RIEMANN_PROBLEM_H

#include <iostream>
#include <cassert>
#include <cmath>

const double GAMMA = 1.4;
//const double K = 1.38064852e-23;
const double K = 1.0;


struct STATE
{
    double u;       //velocity
    double rho;     //density
    double pres;    //pressure

    void print() const
    {
        printf("(%g, %g, %g)\n",u,rho,pres);
    }

    friend std::ostream& operator << (std::ostream& out, const STATE& s)
    {
        s.print();
    }
};

enum class PISTONDIR {LEFT,RIGHT};

double near_piston_soundspeed(double u1, double rho0, double u0, PISTONDIR dir);
//double near_piston_soundspeed(double k, double gamma, double u1, double rho0, double u0, PISTONDIR dir);

double compute_density(double a);
//double compute_density(double k, double gamma, double a);

//These two pressure functions return identical values.
double compute_pressure(double rho0, double p0, double rho1);
//double compute_pressure(double gamma, double rho0, double p0, double rho1);

double compute_pressure(double rho);
//double compute_pressure(double k, double gamma, double rho);

// valid for  (a0 - 0.5T(gamma+1.0)*U)t <= x <= a0*t (for left withdrawal problem)
double compute_fan_velocity(double a0, double u0, double x, double t, PISTONDIR dir);
//double compute_fan_velocity(double gamma, double a0, double u0, double x, double t);

double compute_fan_soundspeed(double a0, double u0, double x, double t, PISTONDIR dir);
//double compute_fan_soundspeed(double gamma, double a0, double u0, double x, double t);




#endif
