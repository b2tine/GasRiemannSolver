#include "riemann_problem.h"

//TODO: get GAMMA from STATE or as function parameter


STATE::STATE(double RHO, double U, double P)
    : rho{RHO}, u{U}, p{P}
{
    if (rho > 0.0)
        a = constant_state_soundspeed(rho,p);
}

STATE::STATE(double RHO, double U, double P, const std::string& ID)
    : STATE{RHO,U,P}
{
    id = ID;
}

STATE::STATE(double RHO, double U, double P, double A)
    : rho{RHO}, u{U}, p{P}, a{A}
{}

STATE::STATE(double RHO, double U, double P, double A, const std::string& ID)
    : STATE{RHO,U,P,A}
{
    id = ID;
}

std::string STATE::printinfo() const
{
    char ostring[250];
    sprintf(ostring,"%5s (%g, %g, %g, %g)",
            ((id.empty()) ? "" : id + " :").c_str(),rho,u,p,a);
    return std::string(ostring);
}


STATE RiemannProblem::operator()(double ksi)
{
    //Locate solution region of ksi = x/t and compute the solution
    if (ksi < slip_slope)
    {
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
                DIRECTION dir = DIRECTION::RIGHT;
                double u_fan = rarefaction_velocity(ksi,dir,sl->u,sl->a);
                double a_fan = rarefaction_soundspeed(ksi,dir,sl->u,sl->a);
                double rho_fan = isentropic_relation_density(a_fan,sl->rho,sl->p);
                double p_fan = isentropic_relation_pressure(a_fan,rho_fan);
                return STATE(rho_fan,u_fan,p_fan,a_fan);
            }
            else
                return *sl_c;
        }
    }
    else
    {
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
                DIRECTION dir = DIRECTION::LEFT;
                double u_fan = rarefaction_velocity(ksi,dir,sr->u,sr->a);
                double a_fan = rarefaction_soundspeed(ksi,dir,sr->u,sr->a);
                double rho_fan = isentropic_relation_density(a_fan,sr->rho,sr->p);
                double p_fan = isentropic_relation_pressure(a_fan,rho_fan);
                return STATE(rho_fan,u_fan,p_fan,a_fan);
            }
            else
                return *sr;
        }
    }
}

void RiemannProblem::solve()
{
    //Solve F(Pslip) = ul_star(P_slip) - ur-star(P_slip) = 0
    Pslip = secantMethod(rpfunc,sl->p,sr->p);
    detectVacuumState();

    //Compute defining characteristics of Riemann Solution
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
        
double LeftCenteredWave(double Pslip, STATE* sl, STATE* sl_center)
{
    if (Pslip > sl->p)
    {
        //LCW is Left Facing Shock Wave (LFS)
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
        //LCW is a GAMMA PLUS Simple Wave (S+)
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
        //RCW is a GAMMA MINUS Simple Wave (S-)
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
        //RCW is a Right Facing Shock Wave (RFS)
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

//TODO: implement variable gamma

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
    /*return (u0*(gamma1-1.0)
            + 2.0*(ksi - sign*a0*(gamma1-1.0)/(gamma0-1.0)))/(gamma1+1.0);*/
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
    /*double a_fan = a0*(gamma1-1.0)/(gamma0-1.0) 
        + sign*(ksi - sign*a0*(gamma1-1.0)/(gamma0-1.0) - u0)*(gamma1-1.0)/(gamma1+1.0);*/
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


