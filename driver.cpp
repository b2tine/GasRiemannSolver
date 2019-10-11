#include "riemann_problem.h"

int randsign()
{
    static bool first = true;
    if (first)
    {
        srand(time(nullptr));
        first = false;
    }
    return 2*(rand() % 2) - 1;
}

double LeftCenteredWave(double Pslip, const STATE& sl)
{
    if (Pslip <= sl.p)
    {
        //TODO: ul_star is rarefaction wave
        return {};
    }
    else
    {
        //ul_star is Left Facing Shock wave
        double ua = sl.u;
        double rhoa = sl.rho;
        double pa = sl.p;
        double pb = Pslip;

        double taua = 1.0/rhoa;
        double taub = behind_state_specific_volume(rhoa,pa,pb);
        
        //M > 0 this case
        double M = std::sqrt((pb - pa)/(taua - taub));
        double ub = ua - (pb - pa)/M;
        return ub;
    }
}

double RightCenteredWave(double Pslip, const STATE& sr)
{
    if (Pslip <= sr.p)
    {
        //TODO: ur_star is rarefaction wave
        return {};
    }
    else
    {
        //ur_star is a Right Facing Shock wave
        double ua = sl.u;
        double rhoa = sl.rho;
        double pa = sl.p;
        double pb = Pslip;

        double taua = 1.0/rhoa;
        double taub = behind_state_specific_volume(rhoa,pa,pb);
        
        //M < 0 this case
        double M = -std::sqrt((pb - pa)/(taua - taub));
        double ub = ua - (pb - pa)/M;
        return ub;
    }
}

struct RP_Function
{
    STATE sl, sr;

    RP_Function(const STATE& sleft, const STATE& sright)
        : sl{sleft}, sr{sright}
    {}

    double operator () (double P) const
    {
        return LeftCenteredWave(P,sl) - RightCenteredWave(P,sr);
    }
};


int main(int argc, char* argv[])
{
    //Given left and right states sl and sr
    
    double ul = 0.0;
    double rhol = 1.225;
    double pl = 100000;

    double ur = -1.0;
    double rhor = 1.1;
    double pr = 85000;

    STATE sl = {ul, rhol, pl};
    STATE sr = {ur, rhor, pr};

    //Since velocities are equal on each side of the contact discontinuity, set
    //    F(P_slip) = ul_star(P_slip) - ur-star(P_slip) = 0
    
    //Need two initial guesses for P_slip to use secant method
    double half_prange = 0.5*std::abs(pl-pr);
    double pg0 = 0.5*(pl+pr);
    double pg1 = pg0 + randsign()*half_prange;

    //Solve F(Pslip) = 0 using the secant method
    double Pslip = secantMethod(RP_Function,pg0,pg1);
    
    
    //TODO: Write functions for ul_star and ur_star, or
    //      LeftCentereredWave() and RightCenteredWave().



    //.... to find the pressure at the
    //contact discontinuity --> use to compute velocities and densities.
    //
    //  For each iteration, compare P to the ahead and behind state pressures to
    //  determine if the state if the LCW/RCW is a shock or rarefaction, and compute
    //  the appropriate value for each.
    //


    return 0;
}
