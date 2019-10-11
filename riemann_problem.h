#ifndef RIEMANN_PROBLEM_H
#define RIEMANN_PROBLEM_H

#include <iostream>
#include <cassert>
#include <cmath>

const double GAMMA = 1.4;


struct STATE
{
    double u;       //velocity
    double rho;     //density
    double p;    //pressure

    void print() const
    {
        printf("(%g, %g, %g)\n",u,rho,p);
    }

    friend std::ostream& operator << (std::ostream& out, const STATE& s)
    {
        s.print();
    }
};


enum class PISTONDIR {LEFT,RIGHT};


double constant_state_soundspeed(double rho, double pres);

double near_piston_soundspeed(double u1, PISTONDIR dir,
                              double u0, double rho0, double pres0);

double near_piston_density(double a1, double rho0, double pres0);

double near_piston_pressure(double a1, double rho1);

double compute_fan_velocity(double x, double t, PISTONDIR dir,
                            double u0, double a0);

double compute_fan_soundspeed(double x, double t, PISTONDIR dir,
                              double u0, double a0);




#endif
