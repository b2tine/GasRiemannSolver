#include "sw_riemann_problem.h"

//TODO: get GAMMA from STATE or as function parameter


STATE::STATE(double U, double H)
    : u{U}, h{H}
{
    computePressure();
    computeCelerity();
}

STATE::STATE(double U, double H, const std::string& ID)
    : STATE{U,H}
{
    id = ID;
}

void STATE::computePressure()
{
    assert(h >= 0.0);
    p = 0.5*G*h*h;
}

void STATE::computeCelerity()
{
    assert(h >= 0.0);
    c = sqrt(G*h);
}

std::string STATE::printinfo() const
{
    char ostring[250];
    sprintf(ostring,"%5s (%g, %g, %g, %g)",
            ((id.empty()) ? "" : id + " :").c_str(),u,h,p,c);
    return std::string(ostring);
}

//Locate solution region of ksi = x/t and compute the solution
STATE RiemannProblem::operator()(double ksi)
{
    if (LCW == WAVETYPE::SHOCK)
    {
        if (ksi < left_shockspeed)
            return *sl;
        else
        {
            if (RCW == WAVETYPE::SHOCK)
            {
                if (ksi < right_shockspeed)
                    return *sl_c;
                else
                    return *sr;
            }
            else if (RCW == WAVETYPE::SIMPLE)
            {
                if (ksi < right_leading_fan_slope)
                    return *sl_c;
                else if (ksi < right_trailing_fan_slope)
                {
                    double u_fan = (sr->u - 2.0*sqrt(G*sr->h) + 2.0*ksi)/3.0;
                    double h_fan = (sr->u - 2.0*sqrt(G*sr->h) - ksi)/3.0;
                    h_fan *= h_fan/G;
                    return STATE(u_fan,h_fan);
                }
                else
                    return *sr;
            }
        }
    }
    else if(LCW == WAVETYPE::SIMPLE)
    {
        if (ksi < left_trailing_fan_slope)
            return *sl;
        else if (ksi < left_leading_fan_slope)
        {
            double u_fan = (sl->u + 2.0*sqrt(G*sl->h) + 2.0*ksi)/3.0;
            double h_fan = (sl->u + 2.0*sqrt(G*sl->h) - ksi)/3.0;
            h_fan *= h_fan/G;
            return STATE(u_fan,h_fan);
        }
        else
        {
            if (RCW == WAVETYPE::SHOCK)
            {
                if (ksi < right_shockspeed)
                    return *sl_c;
                else
                    return *sr;
            }
            else if (RCW == WAVETYPE::SIMPLE)
            {
                if (ksi < right_leading_fan_slope)
                    return *sl_c;
                else if (ksi < right_trailing_fan_slope)
                {
                    double u_fan = (sr->u - 2.0*sqrt(G*sr->h) + 2.0*ksi)/3.0;
                    double h_fan = (sr->u - 2.0*sqrt(G*sr->h) - ksi)/3.0;
                    h_fan *= h_fan/G;
                    return STATE(u_fan,h_fan);
                }
                else
                    return *sr;
            }
        }
    }
}

//TODO: if dry bed dambreak, then only 2 constant states
//      seperated by rarefaction wave. (skip secant method)
void RiemannProblem::solve()
{
    //Solve F(H_ctr) = ul_center(H_ctr) - ur_center(H_ctr) = 0
    H_ctr = secantMethod(rpfunc,sl->h,sr->h);
    if (fabs(sl_c->u) < 1.0e-9)
    {
        sl_c->u = 0.0;
        sr_c->u = 0.0;
    }

    //Compute defining characteristics of Riemann Solution
    if (H_ctr > sl->h)
    {
        //LEFT_FACING_SHOCK
        LCW = WAVETYPE::SHOCK;
        left_shockspeed =
            sl->u - sqrt(0.5*G*sl_c->h*(sl->h + sl_c->h)/sl->h);
    }
    else
    {
        //GAMMA_PLUS_WAVE
        LCW = WAVETYPE::SIMPLE;
        left_trailing_fan_slope = sl->u - sqrt(G*sl->h);
        left_leading_fan_slope = sl_c->u - sqrt(G*sl_c->h);
    }

    if (H_ctr > sr->h)
    {
        //RIGHT_FACING_SHOCK
        RCW = WAVETYPE::SHOCK;
        right_shockspeed =
            sr->u + sqrt(0.5*G*sr_c->h*(sr->h + sr_c->h)/sr->h);
    }
    else
    {
        //GAMMA_MINUS_WAVE
        RCW = WAVETYPE::SIMPLE;
        right_leading_fan_slope = sr_c->u + sqrt(G*sr_c->h);
        right_trailing_fan_slope = sr->u + sqrt(G*sr->h);
    }
}

void RiemannProblem::printStates()
{
    std::cout << "\n";
    std::cout << *sl << "\n";
    std::cout << *sl_c << "\n";
    std::cout << *sr_c << "\n";
    std::cout << *sr << "\n";
}

void RiemannProblem::printWaves()
{
    printf("\n");
    if (LCW == WAVETYPE::SHOCK)
    {
        printf("LCW is a Left Facing Shock:\n");
        printf("\tleft_shockspeed = %g\n",left_shockspeed);
    }
    else if (LCW == WAVETYPE::SIMPLE)
    {
        printf("LCW is a Gamma+ Simple Wave:\n");
        printf("\tleft_trailing_fan_slope = %g\n",left_trailing_fan_slope);
        printf("\tleft_leading_fan_slope = %g\n",left_leading_fan_slope);
    }
    
    printf("\n");
    if (RCW == WAVETYPE::SHOCK)
    {
        printf("RCW is a Right Facing Shock:\n");
        printf("right_shockspeed = %g\n",right_shockspeed);
    }
    else if (RCW == WAVETYPE::SIMPLE)
    {
        printf("RCW is a Gamma- Simple Wave:\n");
        printf("right_leading_fan_slope = %g\n",right_leading_fan_slope);
        printf("right_trailing_fan_slope = %g\n",right_trailing_fan_slope);
    }
}

double LeftCenteredWave(double H, STATE* sl, STATE* sl_center)
{
    double ul, hl;
    double u_lc, h_lc;

    if (H > sl->h)
    {
        //LCW is Left Facing Shock Wave (LFS)
        ul = sl->u;
        hl = sl->h;

        h_lc = H;
        u_lc = ul - (h_lc - hl)*sqrt(0.5*G*(hl + h_lc)/(hl*h_lc));
        
    }
    else
    {
        //LCW is a GAMMA PLUS Simple Wave (S+)
        ul = sl->u;
        hl = sl->h;

        h_lc = H;
        u_lc = ul - 2.0*(sqrt(G*h_lc) - sqrt(G*hl));
 
    }

    //save center state variables
    sl_center->u = u_lc;
    sl_center->h = h_lc;
    sl_center->computePressure();
    sl_center->computeCelerity();
        
    return u_lc;
}

double RightCenteredWave(double H, STATE* sr_center, STATE* sr)
{
    double ur, hr;
    double u_rc, h_rc;

    if (H > sr->h)
    {
        //RCW is a Right Facing Shock Wave (RFS)
        ur = sr->u;
        hr = sr->h;

        hr_c = H;
        u_rc = ur + (h_rc - hr)*sqrt(0.5*G*(hr + h_rc)/(hr*h_rc));
    }
    else
    {
        //RCW is a GAMMA MINUS Simple Wave (S-)
        ur = sr->u;
        hr = sr->h;

        h_rc = H;
        u_rc = ur + 2.0*(sqrt(G*h_rc) - sqrt(G*hr));
    }
    
    //save center state variables
    sr_center->u = u_rc;
    sr_center->h = h_rc;
    sr_center->computePressure();
    sr_center->computeCelerity();
        
    return u_rc;
}
