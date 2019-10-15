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
        //LCW is Left Facing Shock (LFS) Wave
        double ua = sl->u;
        double rhoa = sl->rho;
        double pa = sl->p;
        double pb = Pslip;

        double taua = 1.0/rhoa;
        double taub = behind_state_specific_volume(rhoa,pa,pb);
        double rhob = 1.0/taub;
        
        //M > 0 this case
        double M = std::sqrt((pb - pa)/(taua - taub));
        double ub = ua - (pb - pa)/M;
        
        //save center state variables
        sl_center->u = ub;
        sl_center->rho = rhob;
        sl_center->p = pb;
        sl_center->a = constant_state_soundspeed(rhob,pb);

        return ub;
    }
    else
    {
        //LCW is a GAMMA PLUS (S+) Simple Wave
        double u0 = sl->u;
        double rho0 = sl->rho;
        double p0 = sl->p;
        double a0 = sl->a;

        //near slip line state variables
        double p1 = Pslip;
        double rho1 = rho0*pow(p1/p0,1.0/GAMMA);
        double a1 = constant_state_soundspeed(rho1,p1);
        double u1 = u0 + 2.0*(a0 - a1)/(GAMMA-1.0);

        //save center state variables
        sl_center->u = u1;
        sl_center->rho = rho1;
        sl_center->p = p1;
        sl_center->a = a1;

        return u1;
    }
}

double RightCenteredWave(double Pslip, STATE* sr_center, STATE* sr)
{
    if (Pslip <= sr->p)
    {
        //RCW is a GAMMA MINUS (S-) Simple Wave
        double u0 = sr->u;
        double rho0 = sr->rho;
        double p0 = sr->p;
        double a0 = sr->a;

        //near slip line state variables
        double p1 = Pslip;
        double rho1 = rho0*pow(p1/p0,1.0/GAMMA);
        double a1 = constant_state_soundspeed(rho1,p1);
        double u1 = u0 - 2.0*(a0 - a1)/(GAMMA-1.0);
        
        //save center state variables
        sr_center->u = u1;
        sr_center->rho = rho1;
        sr_center->p = p1;
        sr_center->a = a1;

        return u1;
    }
    else
    {
        //RCW is a Right Facing Shock (RFS) Wave
        double ua = sr->u;
        double rhoa = sr->rho;
        double pa = sr->p;
        double pb = Pslip;

        double taua = 1.0/rhoa;
        double taub = behind_state_specific_volume(rhoa,pa,pb);
        double rhob = 1.0/taub;
        
        //M < 0 this case
        double M = -std::sqrt((pb - pa)/(taua - taub));
        double ub = ua - (pb - pa)/M;
        
        //save center state variables
        sr_center->u = ub;
        sr_center->rho = rhob;
        sr_center->p = pb;
        sr_center->a = constant_state_soundspeed(rhob,pb);

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
    double al = constant_state_soundspeed(rhol,pl);

    double ur = -1.0;
    double rhor = 1.1;
    double pr = 85000;
    double ar = constant_state_soundspeed(rhor,pr);

    STATE sl = {ul, rhol, pl, al};
    STATE sl_c = {0, 0, 0, 0};
    STATE sr_c = {0, 0, 0, 0};
    STATE sr = {ur, rhor, pr, ar};

    //Need two initial guesses for P_slip to use secant method.
    //Just using pl and pr for now, seems better than previous method
    
    //TODO: check if below works better now that secant method bug fixed
    
        /*
        double frac_prange = 0.05*std::abs(pl-pr);
        double pg0 = 0.5*(pl+pr);
        double pg1 = pg0 + randsign()*frac_prange;
        printf("frac_prange = %g\n",frac_prange);
        printf("pg0 = %g\npg1 = %g\n\n",pg0,pg1);
        */

    //Solve F(Pslip) = ul_star(P_slip) - ur-star(P_slip) = 0
    //using the secant method to find the pressure at the contact discontinuity
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

    WAVETYPE LCW, RCW;
    double left_shockspeed;
    double left_trailing_fan_slope;
    double left_leading_fan_slope;
    double right_trailing_fan_slope;
    double right_leading_fan_slope;
    double right_shockspeed;

    if (Pslip > sl.p)
    {
        LCW = WAVETYPE::SHOCK;
        left_shockspeed =
            (sl.rho*sl.u - sl_c.rho*sl_c.u)/(sl.rho - sl_c.rho);
    }
    else
    {
        LCW = WAVETYPE::SIMPLE;
        double al = constant_state_soundspeed(sl.rho,sl.p);
        left_trailing_fan_slope = sl.u - al;

        double al_c = constant_state_soundspeed(sl_c.rho,sl_c.p);
        left_leading_fan_slope = sl_c.u - al_c;
    }

    double slip_slope = sl_c.u;

    if (Pslip < sr.p)
    {
        RCW = WAVETYPE::SIMPLE;
        double ar_c = constant_state_soundspeed(sr_c.rho,sr_c.p);
        right_leading_fan_slope = sr_c.u - ar_c;

        double ar = constant_state_soundspeed(sr.rho,sr.p);
        right_trailing_fan_slope = sr.u - ar;
    }
    else
    {
        RCW = WAVETYPE::SHOCK;
        right_shockspeed =
            (sr.rho*sr.u - sr_c.rho*sr_c.u)/(sr.rho - sr_c.rho);
    }


    //TODO: Given any (x,t), compare ksi = x/t to the slopes
    //      computed above to determine the solution
    
    //sample (x,t) point 
    double x = 0.5;
    double t = 0.02;
    double ksi = x/t;
    double u_riemann;

    if (ksi < slip_slope)
    {
        if (LCW == WAVETYPE::SHOCK)
        {
            if (ksi < left_shockspeed)
                u_riemann = sl.u;
            else
                u_riemann = sl_c.u;
        }
        else
        {
            if (ksi < left_trailing_fan_slope)
                u_riemann = sl.u;
            else if (ksi < left_leading_fan_slope)
            {
                u_riemann = ((GAMMA-1.0)*sl.u + 2.0*(ksi + sl.a))/(GAMMA+1.0); 
            }
            else
                u_riemann = sl_c.u;
        }
    }
    else
    {
        if (RCW == WAVETYPE::SHOCK)
        {
            if (ksi < right_shockspeed)
                u_riemann = sr_c.u;
            else
                u_riemann = sr.u;

        }
        else
        {
            if (ksi < right_leading_fan_slope)
                u_riemann = sr_c.u;
            else if (ksi < right_trailing_fan_slope)
            {
                u_riemann = ((GAMMA-1.0)*sr.u + 2.0*(ksi - sr.a))/(GAMMA+1.0); 
            }
            else
                u_riemann = sr.u;

        }
    }

    std::cout << "u_riemann = " << u_riemann << "\n";

    return 0;
}

