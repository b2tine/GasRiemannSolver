#include "riemann_problem.h"
#include "secant_method.h"

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
        //LCW is a S+ Simple Wave
        double u0 = sl.u;
        double rho0 = sl.rho;
        double p0 = sl.p;
        double a0 = constant_state_soundspeed(rho0,p0);

        //near slip line state variables
        double p1 = Pslip;
        double rho1 = rho0*pow(p1/p0,1.0/GAMMA);
        double a1 = constant_state_soundspeed(rho1,p1);
        double u1 = u0 + 2.0*(a0 - a1)/(GAMMA-1.0);
        return u1;
    }
    else
    {
        //LCW is Left Facing Shock Wave
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
        //RCW is a S- Simple Wave
        double u0 = sr.u;
        double rho0 = sr.rho;
        double p0 = sr.p;
        double a0 = constant_state_soundspeed(rho0,p0);

        //near slip line state variables
        double p1 = Pslip;
        double rho1 = rho0*pow(p1/p0,1.0/GAMMA);
        double a1 = constant_state_soundspeed(rho1,p1);
        double u1 = u0 - 2.0*(a0 - a1)/(GAMMA-1.0);
        return u1;
    }
    else
    {
        //RCW is a Right Facing Shock Wave
        double ua = sr.u;
        double rhoa = sr.rho;
        double pa = sr.p;
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

    //Solve F(Pslip) = 0 using the secant method to
    //find the pressure at the contact discontinuity
    RP_Function F(sl,sr);
    double Pslip = secantMethod(F,pg0,pg1);
    
    std::cout << "sl = " << sl;
    std::cout << "sr = " << sr;
    std::cout << "Pslip = " << Pslip << "\n";
    
    //TODO: compute shocklines, fan regions etc.
    //      provide solution for any (x,t) pair;
    //      locate region of point and compute solution.

    return 0;
}
