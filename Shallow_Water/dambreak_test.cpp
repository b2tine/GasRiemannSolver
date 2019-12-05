#include "sw_riemann_problem.h"


void piston_withdraw_point_locate(
        double x, double t, DIRECTION dir,
        double u1, double u0, double rho0, double p0, double a0);


//rho0, p0, and u0 are state values in the constant region
//far from the moving piston head.
//u1 is the velocity with which the piston head is withdrawn,
//and rho1, p1 are near-piston state variables.

//LEFT piston withdraw creates a "GAMMA_MINUS" simple wave (S-)
//RIGHT piston withdraw creates a "GAMMA_PLUS" simple wave (S+)

int main(int argc, char** argv)
{
    double sign = 1.0;
    DIRECTION dir = DIRECTION::LEFT;
    if (argc > 1)
    {
        if (argv[1][0] == 'r' || argv[1][0] == 'R')
        {
            sign = -1.0;
            dir = DIRECTION::RIGHT;
        }
    }

    //GIVEN:
    double H0 = 1.0;
    //double u1 = sign*U;

    printf("%s Dam Height = %g m\n\n",
            (dir == DIRECTION::LEFT) ? "LEFT" : "RIGHT", H0);

    double u0 = 0.0;        // m/s
    double h0 = H0;         // m

    //double rho0 = 1.225;    // kg/m^3
    //double p0 = 100000;     // kg/m/s^2 (pascals)
    
    //END GIVEN
    double a0 = constant_state_soundspeed(rho0,p0);

    printf("u0 = %g\n",u0);
    printf("h0 = %g\n",h0);

    //printf("rho0 = %g\n",rho0);
    //printf("p0 = %g\n",p0);
    printf("a0 = %g\n\n",a0);

    double u1 = 0.0;
    double h1 = 0.0;

    double a1 = near_piston_soundspeed(u1,dir,u0,rho0,p0);
    double rho1 = isentropic_relation_density(a1,rho0,p0);
    double p1 = isentropic_relation_pressure(a1,rho1);

    printf("u1 = %g\n",u1);
    printf("h1 = %g\n",h1);

    //printf("rho1 = %g\n",rho1);
    //printf("p1 = %g\n",p1);
    //printf("a1 = %g\n\n",a1);

    //sample point to test point location
    double x = sign*3.37;
    double t = 0.01;
    
    piston_withdraw_point_locate(x,t,dir,u1,u0,rho0,p0,a0);
    
    return 0;
}


//TODO: pass in STATE* data structures
void piston_withdraw_point_locate(
        double x, double t, DIRECTION dir,
        double u1, double u0, double rho0, double p0, double a0)
{
    double C = x/t;
    printf("C = %g\n",C);
    if (dir == DIRECTION::LEFT)
    {
        double CplusLeft = 0.5*((GAMMA+1.0)*u1 - (GAMMA-1.0)*u0) + a0;
        double CplusRight = u0 + a0;
        printf("CplusLeft = %g\n",CplusLeft);
        printf("CplusRight = %g\n",CplusRight);

        if (CplusLeft <= C && C <= CplusRight)
        {
            printf("(x,t) = (%g, %g) in fan region, S-\n\n",x,t);
            //TODO: Return u_fan, rho_fan, p_fan.
            double u_fan = rarefaction_velocity_xt(x,t,dir,u0,a0);
            double a_fan = rarefaction_soundspeed_xt(x,t,dir,u0,a0);
            double rho_fan = isentropic_relation_density(a_fan,rho0,p0);
            double p_fan = isentropic_relation_pressure(a_fan,rho_fan);
            
            printf("u_fan = %g\n",u_fan);
            printf("rho_fan = %g\n",rho_fan);
            printf("p_fan = %g\n",p_fan);
            printf("a_fan = %g\n",a_fan);
        }
        else if (C < CplusLeft)
        {
            printf("(x,t) = (%g, %g) in near-piston "
                    "constant region, S1\n\n",x,t);
            //TODO: return u1, rho1, p1
        }
        else
        {
            printf("(x,t) = (%g, %g) in far-piston "
                    "constant region, S0\n\n",x,t);
            //TODO: return u0, rho0, p0
        }
    }
    else
    {
        double CminusLeft = u0 - a0;
        double CminusRight = 0.5*((GAMMA+1.0)*u1 - (GAMMA-1.0)*u0) - a0;
        printf("CminusLeft = %g\n",CminusLeft);
        printf("CminusRight = %g\n",CminusRight);
        
        if (CminusLeft <= C && C <= CminusRight)
        {
            printf("(x,t) = (%g, %g) in fan region, S+\n\n",x,t);
            //TODO: Return u_fan, rho_fan, p_fan.
            double u_fan = rarefaction_velocity_xt(x,t,dir,u0,a0);
            double a_fan = rarefaction_soundspeed_xt(x,t,dir,u0,a0);
            double rho_fan = isentropic_relation_density(a_fan,rho0,p0);
            double p_fan = isentropic_relation_pressure(a_fan,rho_fan);
            
            printf("u_fan = %g\n",u_fan);
            printf("rho_fan = %g\n",rho_fan);
            printf("p_fan = %g\n",p_fan);
            printf("a_fan = %g\n",a_fan);
        }
        else if (CminusRight < C)
        {
            printf("(x,t) = (%g, %g) in near-piston "
                    "constant region, S1\n\n",x,t);
            //TODO: return u1, rho1, p1
        }
        else
        {
            printf("(x,t) = (%g, %g) in far-piston "
                    "constant region, S0\n\n",x,t);
            //TODO: return u0, rho0, p0
        }
    }
}
