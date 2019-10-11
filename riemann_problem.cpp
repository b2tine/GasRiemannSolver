#include "riemann_problem.h"


double near_piston_soundspeed(double u1, PISTONDIR dir,
        double u0, double rho0, double pres0)
       
{
    double sign = ((dir == PISTONDIR::LEFT) ? 1.0 : -1.0);
    double a0 = constant_state_soundspeed(rho0,pres0);
    double a = a0 + sign*0.5*(GAMMA-1.0)*(u1-u0);
    
    if (a <= 0.0)
    {
        //TODO: throw an exception
        printf("Vacuum State Detected near piston head\n");
        exit(1);
    }
    return a;
}

double constant_state_soundspeed(double rho, double pres)
{
    return std::sqrt(GAMMA*pres/rho);
}

//a1 is the near piston soundspeed
double near_piston_density(double a1, double rho0, double pres0)
{
    double arg = a1*a1*pow(rho0,GAMMA)/GAMMA/pres0;
    return pow(arg,1.0/(GAMMA-1.0));
}

//a1 is the near piston soundspeed, rho1 is near piston density
double near_piston_pressure(double a1, double rho1)
{
    return a1*a1*rho1/GAMMA;
}

//TODO: make pure functions of x/t (combine x,t into single argument)?
double compute_fan_velocity(double x, double t, PISTONDIR dir,
                            double u0, double a0)
{
    assert(t > 0.0);
    double sign = ((dir == PISTONDIR::LEFT) ? 1.0 : -1.0);
    return (2.0*(x/t - sign*a0) + u0*(GAMMA-1.0))/(GAMMA+1.0);
}

double compute_fan_soundspeed(double x, double t, PISTONDIR dir,
                              double u0, double a0)
{
    assert(t > 0.0);
    double sign = ((dir == PISTONDIR::LEFT) ? 1.0 : -1.0);
    double a_fan = a0 + sign*(x/t - sign*a0 - u0)*(GAMMA-1.0)/(GAMMA+1.0);
    if (a_fan <= 0.0)
    {
        //TODO: throw an exception
        printf("Vacuum State Detected in fan region\n");
        exit(1);
    }
    return a_fan;
}

//TODO: Write a Left/Right CenteredWave() function per assignment instructions.
//      Should just call the fan solutions for corresponding directions? 


