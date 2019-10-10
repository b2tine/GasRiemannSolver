#include "riemann_problem.h"



int main(int argc, char* argv[])
{
    //double k = 1.0;
    //double gamma = 1.4;

    //Given behind and ahead states (sb and sa)
    
    double u_b = 0.0;
    double rho_b = 1.225;
    double pres_b = compute_pressure(rho_b);

    double u_a = -1.0;
    double rho_a = 1.1;
    double pres_a = compute_pressure(rho_a);

    STATE sb = {u_b, rho_b, pres_b};
    STATE sa = {u_a, rho_a, pres_a};

    //Since velocities are equal on each side of the contact discontinuity, set
    //    F(P) = ul_star(P) - ur-star(P)
    
    //TODO: Write functions for ul_star and ur_star, or
    //      LeftCentereredWave() and RightCenteredWave().
    //
    //Solve F(P) = 0 using the secant method to find the pressure at the
    //contact discontinuity --> use to compute velocities and densities.
    //
    //  For each iteration, compare P to the ahead and behind state pressures to
    //  determine if the state if the LCW/RCW is a shock or rarefaction, and compute
    //  the appropriate value for each.
    //


    return 0;
}
