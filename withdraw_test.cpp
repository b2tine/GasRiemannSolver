#include "riemann_problem.h"

//rho0, pres0, and u0 are state values in the constant region
//far from the moving piston head.
//U is the velocity with which the piston head is withdrawn.

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

    //double k = 1.0;
    //double gamma = 1.4;

    double rho0 = 1.225;
    double u0 = 0.0;
    
    double Uleft = -1.0;
    //double Uleft = -0.000000000001;
    double U = sign*Uleft;
    printf("U = %g\n",U);

    double pres0 = compute_pressure(rho0);
    double a0 = sqrt(K*GAMMA*pow(rho0,GAMMA-1.0));

    printf("a0 = %g\n",a0);
    printf("rho0 = %g\n",rho0);
    printf("pres0 = %g\n\n",pres0);

    double a1 = near_piston_soundspeed(U,rho0,u0,dir);
    double rho1 = compute_density(a1);
    double pres1 = compute_pressure(rho1);

    printf("a1 = %g\n",a1);
    printf("rho1 = %g\n",rho1);
    printf("pres1 = %g\n\n",pres1);

    
    //sample point in fan region
    double x = -sign*0.1;
    double t = 0.1;
    
    double u_fan = compute_fan_velocity(a0,u0,x,t,dir);
    double a_fan = compute_fan_soundspeed(a0,u0,x,t,dir);
    
    printf("u_fan = %g\n",u_fan);
    printf("a_fan = %g\n",a_fan);

    return 0;
}


