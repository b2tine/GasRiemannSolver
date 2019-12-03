#include "riemann_problem.h"

//TODO: get GAMMA from STATE or as function parameter


STATE::STATE(double U, double H)
    : u{U}, h{H}
{}

STATE::STATE(double U, double H, const std::string& ID)
    : STATE{U,H}
{
    id = ID;
}

/*
STATE::STATE(double U, double H, double P, double A)
    : u{U}, h{H}, p{P}, a{A}
{}

STATE::STATE(double U, double H, double P, double A, const std::string& ID)
    : STATE{U,H,P,A}
{
    id = ID;
}
*/

/*
void STATE::computeSoundSpeed()
{
    a = constant_state_soundspeed(rho,p);
}
*/

std::string STATE::printinfo() const
{
    char ostring[250];
    //sprintf(ostring,"%5s (%g, %g, %g, %g)",
      //      ((id.empty()) ? "" : id + " :").c_str(),u,h,p,a);

    sprintf(ostring,"%5s (%g, %g)",
            ((id.empty()) ? "" : id + " :").c_str(),u,h);
    return std::string(ostring);
}


STATE RiemannProblem::operator()(double ksi)
{
    //Locate solution region of ksi = x/t and compute the solution
    if (LCW == WAVETYPE::SHOCK)
    {
        if (ksi < left_shockspeed)
            return *sl;
        else
            return *sl_c;
    }
    else
    {
        if (ksi < left_trailing_fan_slope)
            return *sl;
        else if (ksi < left_leading_fan_slope)
        {
            double u_fan = (sl->u + 2.0*sqrt(G*sl->h) + 2.0*ksi)/3.0;
            double h_fan = (2.0*sqrt(G*sl->h) + sl->u - ksi)/3.0;
            h_fan *= h_fan/G;
            return STATE(u_fan,h_fan);
        }
        else
            return *sl_c;
    }
    
    if (RCW == WAVETYPE::SHOCK)
    {
        if (ksi < right_shockspeed)
            return *sr_c;
        else
            return *sr;
    }
    else
    {
        if (ksi < right_leading_fan_slope)
            return *sr_c;
        else if (ksi < right_trailing_fan_slope)
        {
            double u_fan = (sr->u - 2.0*sqrt(G*sr->h) + 2.0*ksi)/3.0;
            double h_fan = (2.0*sqrt(G*sr->h) - sr->u + ksi)/3.0;
            h_fan *= h_fan/G;
            return STATE(u_fan,h_fan);
        }
        else
            return *sr;
    }
}

void RiemannProblem::solve()
{
    //sl->computeSoundSpeed();
    //sr->computeSoundSpeed();

    //Solve F(H_ctr) = ul_center(H_ctr) - ur_center(H_ctr) = 0
    H_ctr= secantMethod(rpfunc,sl->h,sr->h);
    //detectVacuumState();

    //Compute defining characteristics of Riemann Solution
    if (H_ctr > sl->h)
    {
        //LEFT_FACING_SHOCK
        LCW = WAVETYPE::SHOCK;
        left_shockspeed =
            sl->u - sqrt(0.5*G*H_ctr*(sl->h + H_ctr)/sl->h);
    }
    else
    {
        //GAMMA_PLUS_WAVE
        LCW = WAVETYPE::SIMPLE;
        left_trailing_fan_slope = sl->u - sqrt(G*sl->h);
        left_leading_fan_slope = sl_c->u - sqrt(G*H_ctr);
    }

    if (H_ctr < sr->h)
    {
        //GAMMA_MINUS_WAVE
        RCW = WAVETYPE::SIMPLE;
        right_leading_fan_slope = sr_c->u + sqrt(G*H_ctr);
        right_trailing_fan_slope = sr->u + sqrt(G*sr->h);
    }
    else
    {
        //RIGHT_FACING_SHOCK
        RCW = WAVETYPE::SHOCK;
        right_shockspeed =
            sr->u + sqrt(0.5*G*H_ctr*(H_ctr + sr->h)/sr->h);
    }
}

void RiemannProblem::printStates()
{
    std::cout << *sl << "\n";
    std::cout << *sl_c << "\n";
    std::cout << *sr_c << "\n";
    std::cout << *sr << "\n";
}

/*
void RiemannProblem::detectVacuumState()
{
    std::vector<STATE*> vacstates;
    for (STATE* s : {sl, sl_c, sr_c, sr})
    {
        if (s->a <= 0.0)
            vacstates.push_back(s);
    }

    if (!vacstates.empty())
        throw VacuumStateException("",vacstates);
}
*/
        
double LeftCenteredWave(double Hctr, STATE* sl, STATE* sl_center)
{
    if (Hctr > sl->h)
    {
        //LCW is Left Facing Shock Wave (LFS)
        double ua = sl->u;
        double ha = sl->h;

        double hb = Hctr;
        double ub = ua + (ha - hb)*sqrt(0.5*G*(ha + hb)/(ha*hb));
        
        //save center state variables
        sl_center->u = ub;
        sl_center->h = hb;
        
        return ub;
    }
    else
    {
        //LCW is a GAMMA PLUS Simple Wave (S+)
        double ul = sl->u;
        double hl = sl->h;

        double hl_c = Hctr;
        double ul_c = ul + 2.0*(sqrt(G*hl) - sqrt(G*hl_c));

        //save center state variables
        sl_center->u = ul_c;
        sl_center->h = hl_c;
        
        return ul_c;
    }
}

double RightCenteredWave(double Hctr, STATE* sr_center, STATE* sr)
{
    if (Hctr <= sr->h)
    {
        //RCW is a GAMMA MINUS Simple Wave (S-)
        double ur = sr->u;
        double hr = sr->h;

        double hr_c = Hctr;
        double ur_c = ur + 2.0*(sqrt(G*hr_c) - sqrt(G*hr));

        //save center state variables
        sr_center->u = ur_c;
        sr_center->h = hr_c;
        
        return ur_c;
    }
    else
    {
        //RCW is a Right Facing Shock Wave (RFS)
        double ua = sr->u;
        double ha = sr->h;

        double hb = Hctr;
        double ub = ua + (ha - hb)*sqrt(0.5*G*(ha + hb)/(ha*hb));
        
        //save center state variables
        sr_center->u = ub;
        sr_center->h = hb;
        
        return ub;
    }
}

/*

//SHOCK WAVE FUNCTIONS

double behind_state_specific_volume(double rhoa, double pa, double pb)
{
    double GP = GAMMA + 1.0;
    double GM = GAMMA - 1.0;
    return (GP*pa + GM*pb)/(GP*pb + GM*pa)/rhoa;
    //TODO: from Hugoniot function calculate for variable gamma
}

double behind_state_pressure(double rhoa, double pa, double rhob)
{
    double GP = GAMMA + 1.0;
    double GM = GAMMA - 1.0;
    return (GP/rhoa - GM/rhob)/(GP/rhob - GM/rhoa)*pa;
    //TODO: from Hugoniot function calculate for variable gamma
}


//SIMPLE WAVE FUNCTIONS

double constant_state_soundspeed(double rho, double pres)
{
    return std::sqrt(GAMMA*pres/rho);
}

double near_piston_soundspeed(double u1, DIRECTION dir,
        double u0, double rho0, double pres0)
       
{
    double sign = ((dir == DIRECTION::LEFT) ? 1.0 : -1.0);
    double a0 = constant_state_soundspeed(rho0,pres0);
    double a = a0 + sign*0.5*(GAMMA-1.0)*(u1-u0);
    //Variable Gamma:
    //double a = (a0/(gamma0-1.0) + sign*0.5*(u1-u0))/(gamma1-1.0);
    if (a <= 0.0)
        throw VacuumStateException("near piston head");
    return a;
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

double rarefaction_velocity(double ksi, DIRECTION dir, double u0, double a0)
{
    double sign = ((dir == DIRECTION::LEFT) ? 1.0 : -1.0);
    return (u0*(GAMMA-1.0) + 2.0*(ksi - sign*a0))/(GAMMA+1.0);
    //Variable Gamma:
    //return (u0*(gamma1-1.0)
      //      + 2.0*(ksi - sign*a0*(gamma1-1.0)/(gamma0-1.0)))/(gamma1+1.0);
}

double rarefaction_velocity_xt(double x, double t,
                               DIRECTION dir, double u0, double a0)
{
    assert(t > 0.0);
    return rarefaction_velocity(x/t,dir,u0,a0);
}

double rarefaction_soundspeed(double ksi, DIRECTION dir, double u0, double a0)
{
    double sign = ((dir == DIRECTION::LEFT) ? 1.0 : -1.0);
    double a_fan = a0 + sign*(ksi - sign*a0 - u0)*(GAMMA-1.0)/(GAMMA+1.0);
    //Variable Gamma:
    //double a_fan = a0*(gamma1-1.0)/(gamma0-1.0) 
     //   + sign*(ksi - sign*a0*(gamma1-1.0)/(gamma0-1.0) - u0)*(gamma1-1.0)/(gamma1+1.0);
    if (a_fan <= 0.0)
        throw VacuumStateException("in fan region");
    return a_fan;
}

double rarefaction_soundspeed_xt(double x, double t,
                                 DIRECTION dir, double u0, double a0)
{
    assert(t > 0.0);
    return rarefaction_soundspeed(x/t,dir,u0,a0);
}
*/

