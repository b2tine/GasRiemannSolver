#include "riemann_problem.h"


double compute_ul_star(double Pslip, const STATE& sl)
{
    if (Pslip < sl.p)
    {
        //TODO: ur_star is rarefaction solution
    }
    else
    {
        //TODO: ur_star is shock solution
    }

    return {};

}

double compute_ur_star(double Pslip, const STATE& sr)
{
    if (Pslip < sr.p)
    {
        //TODO: ur_star is rarefaction solution
    }
    else
    {
        //TODO: ur_star is shock solution
    }

    return {};
}

int main(int argc, char* argv[])
{
    //double k = 1.0;
    //double gamma = 1.4;

    //Given behind and ahead states (sb and sa)
    //  use left and right for time being .. sl and sr
    
    double u_l = 0.0;
    double rho_l = 1.225;
    double pres_l = compute_pressure(rho_l);

    double u_r = -1.0;
    double rho_r = 1.1;
    double pres_r = compute_pressure(rho_r);

    STATE sl = {u_l, rho_l, pres_l};
    STATE sr = {u_r, rho_r, pres_r};

    //Since velocities are equal on each side of the contact discontinuity, set
    //    F(P_slip) = ul_star(P_slip) - ur-star(P_slip)
    
    //Need two initial guesses for P_slip
    
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