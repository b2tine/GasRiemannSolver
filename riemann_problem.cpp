#include "riemann_problem.h"


double near_piston_soundspeed(double u1, double rho0, double u0, PISTONDIR dir)
{
    double sign = ((dir == PISTONDIR::LEFT) ? 1.0 : -1.0);
    double a0 = sqrt(K*GAMMA*pow(rho0,GAMMA-1.0));
    double a = a0 + sign*0.5*(GAMMA-1.0)*(u1-u0);
    
    if (a <= 0.0)
    {
        //TODO: how should this be handled, throw exception/error?
        printf("Vacuum State Detected near piston head\n");
        exit(1);
    }
    return a;
}

double compute_density(double a)
{
    return pow(a*a/K/GAMMA,1.0/(GAMMA-1.0));
}

//These two pressure functions return identical values.
double compute_pressure(double rho0, double p0, double rho1)
{
    return p0*pow(rho1/rho0,GAMMA);
}

double compute_pressure(double rho)
{
    return K*pow(rho,GAMMA);
}


//TODO: make pure functions of x/t (combine x,t into single argument)
double compute_fan_velocity(double a0, double u0, double x, double t, PISTONDIR dir)
{
    assert(t > 0.0);
    double sign = ((dir == PISTONDIR::LEFT) ? 1.0 : -1.0);
    return (2.0*(x/t - sign*a0) + u0*(GAMMA-1.0))/(GAMMA+1.0);
    //return (2.0*(x/t - a0) + u0*(GAMMA-1.0))/(GAMMA+1.0);
}

double compute_fan_soundspeed(double a0, double u0, double x, double t, PISTONDIR dir)
{
    assert(t > 0.0);
    double sign = ((dir == PISTONDIR::LEFT) ? 1.0 : -1.0);
    double a_fan = a0 + sign*(x/t - sign*a0 - u0)*(GAMMA-1.0)/(GAMMA+1.0);
    if (a_fan <= 0.0)
    {
        //TODO: how should this be handled, throw exception/error?
        printf("Vacuum State Detected in fan region\n");
        //exit(1);
    }
    return a_fan;
}



