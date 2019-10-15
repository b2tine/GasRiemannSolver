#ifndef RIEMANN_PROBLEM_H
#define RIEMANN_PROBLEM_H

#include "secant_method.h"

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <ctime>
#include <cmath>

const double GAMMA = 1.4;

enum class WAVETYPE {SHOCK,SIMPLE};

struct STATE
{
    double u;       //velocity
    double rho;     //density
    double p;       //pressure
    double a;       //soundspeed
    //double gamma;

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


//SIMPLE WAVE FUNCTIONS

double constant_state_soundspeed(double rho, double pres);

double near_piston_soundspeed(double u1, PISTONDIR dir,
                              double u0, double rho0, double pres0);

double isentropic_relation_density(double a1, double rho0, double pres0);

double isentropic_relation_pressure(double a1, double rho1);

double rarefaction_velocity(double x, double t, PISTONDIR dir,
                            double u0, double a0);

double rarefaction_soundspeed(double x, double t, PISTONDIR dir,
                              double u0, double a0);

//SHOCK WAVE FUNCTIONS

double behind_state_pressure(double rhoa, double pa, double rhob);

double behind_state_specific_volume(double rhoa, double pa, double pb);



#endif
