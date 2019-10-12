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

double LeftCenteredWave(double Pslip, STATE* sl, STATE* sl_center)
{
    if (Pslip > sl->p)
    {
        //LCW is Left Facing Shock Wave
        double ua = sl->u;
        double rhoa = sl->rho;
        double pa = sl->p;
        double pb = Pslip;

        double taua = 1.0/rhoa;
        double taub = behind_state_specific_volume(rhoa,pa,pb);
        
        //M > 0 this case
        double M = std::sqrt((pb - pa)/(taua - taub));
        double ub = ua - (pb - pa)/M;
        
        //save center state variables
        sl_center->u = ub;
        sl_center->rho = 1.0/taub;
        sl_center->p = pb;

        return ub;
    }
    else
    {
        //LCW is a S+ Simple Wave
        double u0 = sl->u;
        double rho0 = sl->rho;
        double p0 = sl->p;
        double a0 = constant_state_soundspeed(rho0,p0);

        //near slip line state variables
        double p1 = Pslip;
        double rho1 = rho0*pow(p1/p0,1.0/GAMMA);
        double a1 = constant_state_soundspeed(rho1,p1);
        double u1 = u0 + 2.0*(a0 - a1)/(GAMMA-1.0);

        //save center state variables
        sl_center->u = u1;
        sl_center->rho = rho1;
        sl_center->p = p1;

        return u1;
    }
}

double RightCenteredWave(double Pslip, STATE* sr_center, STATE* sr)
{
    if (Pslip <= sr->p)
    {
        //RCW is a S- Simple Wave
        double u0 = sr->u;
        double rho0 = sr->rho;
        double p0 = sr->p;
        double a0 = constant_state_soundspeed(rho0,p0);

        //near slip line state variables
        double p1 = Pslip;
        double rho1 = rho0*pow(p1/p0,1.0/GAMMA);
        double a1 = constant_state_soundspeed(rho1,p1);
        double u1 = u0 - 2.0*(a0 - a1)/(GAMMA-1.0);
        
        //save center state variables
        sr_center->u = u1;
        sr_center->rho = rho1;
        sr_center->p = p1;

        return u1;
    }
    else
    {
        //RCW is a Right Facing Shock Wave
        double ua = sr->u;
        double rhoa = sr->rho;
        double pa = sr->p;
        double pb = Pslip;

        double taua = 1.0/rhoa;
        double taub = behind_state_specific_volume(rhoa,pa,pb);
        
        //M < 0 this case
        double M = -std::sqrt((pb - pa)/(taua - taub));
        double ub = ua - (pb - pa)/M;
        
        //save center state variables
        sr_center->u = ub;
        sr_center->rho = 1.0/taub;
        sr_center->p = pb;

        return ub;
    }
}

struct RP_Function
{
    STATE *sl, *sl_cen, *sr_cen, *sr;

    RP_Function(STATE* sL, STATE* sLC, STATE* sRC, STATE* sR)
        : sl{sL}, sl_cen{sLC}, sr_cen{sRC}, sr{sR}
    {}

    double operator () (double P) const
    {
        return LeftCenteredWave(P,sl,sl_cen)
            - RightCenteredWave(P,sr_cen,sr);
    }
};


int main(int argc, char* argv[])
{
    //Given left and right states sl and sr:
    double ul = 0.0;
    double rhol = 1.225;
    double pl = 100000;

    double ur = -1.0;
    double rhor = 1.1;
    double pr = 85000;

    STATE sl = {ul, rhol, pl};
    STATE sl_c = {0, 0, 0};
    STATE sr_c = {0, 0, 0};
    STATE sr = {ur, rhor, pr};

    //Since velocities are equal on each side of the contact discontinuity, set
    //    F(P_slip) = ul_star(P_slip) - ur-star(P_slip) = 0
    
    //Need two initial guesses for P_slip to use secant method.
    //Just using pl and pr for now, seems better than previous method
    
    //TODO: The order of passing pl and pr as initial guesses matters.
    //      May need to consider velocity and density information in
    //      making this decision. Try using entropy conditions and
    //      properties of the waves
    
        /*
        double frac_prange = 0.05*std::abs(pl-pr);
        double pg0 = 0.5*(pl+pr);
        double pg1 = pg0 + randsign()*frac_prange;
        printf("frac_prange = %g\n",frac_prange);
        printf("pg0 = %g\npg1 = %g\n\n",pg0,pg1);
        */

    //Solve F(Pslip) = 0 using the secant method to
    //find the pressure at the contact discontinuity
    RP_Function F(&sl,&sl_c,&sr_c,&sr);
    double Pslip = secantMethod(F,pl,pr);
    //double Pslip = secantMethod(F,pg0,pg1);
    
    std::cout << "Pslip = " << Pslip << "\n";

    std::cout << "sl = " << sl;
    std::cout << "sl_c = " << sl_c;
    std::cout << "sr_c = " << sr_c;
    std::cout << "sr = " << sr;


    //TODO: make this a function and return a structure
    //      contain the riemann problem solution along with
    //      characteristic slopes seperating them

    if (Pslip > sl.p)
    {
        double left_shockspeed =
            (sl.rho*sl.u - sl_c.rho*sl_c.u)/(sl.rho - sl_c.rho);
    }
    else
    {
        double al = constant_state_soundspeed(sl.rho,sl.p);
        double left_trailing_fan_slope = sl.u - al;

        double al_c = constant_state_soundspeed(sl_c.rho,sl_c.p);
        double left_leading_fan_slope = sl_c.u - al_c;
    }

    double slip_slope = sl_c.u;

    if (Pslip < sr.p)
    {
        double ar_c = constant_state_soundspeed(sr_c.rho,sr_c.p);
        double right_leading_fan_slope = sr_c.u - ar_c;

        double ar = constant_state_soundspeed(sr.rho,sr.p);
        double right_trailing_fan_slope = sr.u - ar;
    }
    else
    {
        double right_shockspeed =
            (sr.rho*sr.u - sr_c.rho*sr_c.u)/(sr.rho - sr_c.rho);
    }

    //sample (x,t) point 
    double x = 0.5;
    double t = 0.02;
    double ksi = x/t;

    //TODO: Given any (x,t), compare ksi = x/t to the slopes
    //      computed above to determine

    //pg 11 of lecture_note4.pdf has comparisons 

    return 0;
}
