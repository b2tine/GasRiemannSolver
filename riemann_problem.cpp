#include "riemann_problem.h"

//TODO: get GAMMA from STATE or as function parameter

void RiemannProblem::solve()
{
    Pslip = secantMethod(rpfunc,sl->p,sr->p);

    //TODO: detect vacuum state
    if (Pslip > sl->p)
    {
        LCW = WAVETYPE::SHOCK;
        left_shockspeed =
            (sl->rho*sl->u - sl_c->rho*sl_c->u)/(sl->rho - sl_c->rho);
    }
    else
    {
        LCW = WAVETYPE::SIMPLE;
        left_trailing_fan_slope = sl->u - sl->a;
        left_leading_fan_slope = sl_c->u - sl_c->a;
    }

    slip_slope = sl_c->u;

    if (Pslip < sr->p)
    {
        RCW = WAVETYPE::SIMPLE;
        right_leading_fan_slope = sr_c->u - sr_c->a;
        right_trailing_fan_slope = sr->u - sr->a;
    }
    else
    {
        RCW = WAVETYPE::SHOCK;
        right_shockspeed =
            (sr->rho*sr->u - sr_c->rho*sr_c->u)/(sr->rho - sr_c->rho);
    }
}
        
double RiemannProblem::operator () (double ksi)
{
    double u_riemann;

    if (ksi < slip_slope)
    {
        if (LCW == WAVETYPE::SHOCK)
        {
            if (ksi < left_shockspeed)
                u_riemann = sl->u;
            else
                u_riemann = sl_c->u;
        }
        else
        {
            if (ksi < left_trailing_fan_slope)
                u_riemann = sl->u;
            else if (ksi < left_leading_fan_slope)
            {
                u_riemann = ((GAMMA-1.0)*sl->u + 2.0*(ksi + sl->a))/(GAMMA+1.0); 
            }
            else
                u_riemann = sl_c->u;
        }
    }
    else
    {
        if (RCW == WAVETYPE::SHOCK)
        {
            if (ksi < right_shockspeed)
                u_riemann = sr_c->u;
            else
                u_riemann = sr->u;

        }
        else
        {
            if (ksi < right_leading_fan_slope)
                u_riemann = sr_c->u;
            else if (ksi < right_trailing_fan_slope)
            {
                u_riemann = ((GAMMA-1.0)*sr->u + 2.0*(ksi - sr->a))/(GAMMA+1.0); 
            }
            else
                u_riemann = sr->u;
        }
    }

    return u_riemann;
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


//SIMPLE WAVE FUNCTIONS

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

//a1 is the near piston soundspeed, or a soundspeed in the rarefaction fan.
double isentropic_relation_density(double a1, double rho0, double pres0)
{
    double arg = a1*a1*pow(rho0,GAMMA)/GAMMA/pres0;
    return pow(arg,1.0/(GAMMA-1.0));
}

//a1 is the near piston soundspeed, rho1 is near piston density.
//Or a1 and rho1 are state values along a characteristic in the rarefaction fan.
double isentropic_relation_pressure(double a1, double rho1)
{
    return a1*a1*rho1/GAMMA;
}

//TODO: make pure functions of x/t (combine x,t into single argument)?
double rarefaction_velocity(double x, double t, PISTONDIR dir,
                            double u0, double a0)
{
    assert(t > 0.0);
    double sign = ((dir == PISTONDIR::LEFT) ? 1.0 : -1.0);
    return (2.0*(x/t - sign*a0) + u0*(GAMMA-1.0))/(GAMMA+1.0);
}

double rarefaction_soundspeed(double x, double t, PISTONDIR dir,
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


//SHOCK WAVE FUNCTIONS

double behind_state_pressure(double rhoa, double pa, double rhob)
{
    double GP = GAMMA + 1.0;
    double GM = GAMMA - 1.0;
    return (GP/rhoa - GM/rhob)/(GP/rhob - GM/rhoa)*pa;
}

double behind_state_specific_volume(double rhoa, double pa, double pb)
{
    double GP = GAMMA + 1.0;
    double GM = GAMMA - 1.0;
    return (GP*pa + GM*pb)/(GP*pb + GM*pa)/rhoa;
}



