#include "riemann_problem.h"


void piston_withdraw_point_locate(double x, double t);


//rho0, p0, and u0 are state values in the constant region
//far from the moving piston head.
//u1 is the velocity with which the piston head is withdrawn.

int main(int argc, char** argv)
{
    double sign = 1.0;
    PISTONDIR dir = PISTONDIR::LEFT;
    if (argc > 1)
    {
        if (argv[1][0] == 'r' || argv[1][0] == 'R')
        {
            sign = -1.0;
            dir = PISTONDIR::RIGHT;
        }
    }

    //double GAMMA = 1.4;
    
    //GIVEN:
    double Uleft = -1.0;    // m/s
    double u1 = sign*Uleft;

    printf("%s piston withdrawal speed U = %g m/s\n\n",
            (dir == PISTONDIR::LEFT) ? "LEFT" : "RIGHT", u1);

    double u0 = 0.0;        // m/s
    double rho0 = 1.225;    // kg/m^3
    double p0 = 100000;     // kg/m/s^2 (pascals)
    double a0 = constant_state_soundspeed(rho0,p0);

    printf("u0 = %g\n",u0);
    printf("rho0 = %g\n",rho0);
    printf("p0 = %g\n",p0);
    printf("a0 = %g\n\n",a0);

    double a1 = near_piston_soundspeed(u1,dir,u0,rho0,p0);
    double rho1 = near_piston_density(a1,rho0,p0);
    double p1 = near_piston_pressure(a1,rho1);

    //sample point in fan region
    double x = -sign*0.1;
    double t = 0.1;
    
    piston_withdraw_point_locate(x,t);
    
    double u_fan = compute_fan_velocity(x,t,dir,u0,a0);
    double a_fan = compute_fan_soundspeed(x,t,dir,u0,a0);
    
    printf("u1 = %g\n",u1);
    printf("rho1 = %g\n",rho1);
    printf("p1 = %g\n",p1);
    printf("a1 = %g\n",a1);

    printf("u_fan = %g\n",u_fan);
    printf("a_fan = %g\n",a_fan);

    return 0;
}


void piston_withdraw_point_locate(double x, double t)
{
    double C = x/t;
    if (dir == PISTONDIR::LEFT)
    {
        double CplusLeft = 0.5*((GAMMA+1.0)*u1 - (GAMMA-1.0)*u0) + a0;
        double CplusRight = u0 + a0;
        if (CplusLeft <= C && C <= CplusRight)
        {
            printf("(x,t) = (%g, %g) in LCW fan region, S-\n\n",x,t);
        }
        else if (C < CplusLeft)
        {
            printf("(x,t) = (%g, %g) in LCW near-piston "
                    "constant region, S1\n\n",x,t);
        }
        else
        {
            printf("(x,t) = (%g, %g) in LCW far-piston "
                    "constant region, S0\n\n",x,t);
        }
    }
    else
    {
        double CminusLeft = u0 - a0;
        double CminusRight = 0.5*((GAMMA+1.0)*u1 - (GAMMA-1.0)*u0) - a0;
        if (CminusLeft <= C && C <= CminusRight)
        {
            printf("(x,t) = (%g, %g) in RCW fan region, S+\n\n",x,t);
        }
        else if (CminusRight < C)
        {
            printf("(x,t) = (%g, %g) in RCW near-piston "
                    "constant region, S1\n\n",x,t);
        }
        else
        {
            printf("(x,t) = (%g, %g) in RCW far-piston "
                    "constant region, S0\n\n",x,t);
        }
    }
}

