#ifndef RIEMANN_PROBLEM_H
#define RIEMANN_PROBLEM_H

#include "secant_method.h"

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <limits>

const double HUGE = std::numeric_limits<double>::max();


const double GAMMA = 1.4;

enum class WAVETYPE {SHOCK,SIMPLE};

struct STATE
{
    double u;       //velocity
    double rho;     //density
    double p;       //pressure
    double a;       //soundspeed
    //double gamma;

    void print() const
    {
        printf("(%g, %g, %g)\n",u,rho,p);
    }

    friend std::ostream& operator << (std::ostream& out, const STATE& s)
    {
        s.print();
    }
};


double LeftCenteredWave(double Pslip, STATE* sl, STATE* sl_center);
double RightCenteredWave(double Pslip, STATE* sr_center, STATE* sr);

/*
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
*/

class RiemannProblem
{
    public:
    
        RiemannProblem(STATE* sL, STATE* sLC, STATE* sRC, STATE* sR)
            : sl{sL}, sl_c{sLC}, sr_c{sRC}, sr{sR}, rpfunc{sL,sLC,sRC,sR}
        {
            this->solve();
        }

        double operator () (double ksi);

        double operator () (double x, double t)
        {
            return this->operator()(x/t);
        }

    private:

        STATE *sl, *sl_c, *sr_c, *sr;
        
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

        } rpfunc;

        double Pslip;
        WAVETYPE LCW, RCW;

        double left_shockspeed {HUGE};
        double left_trailing_fan_slope {HUGE};
        double left_leading_fan_slope {HUGE};
        double slip_slope {HUGE};
        double right_trailing_fan_slope {-HUGE};
        double right_leading_fan_slope {-HUGE};
        double right_shockspeed {-HUGE};

        void solve();
};



enum class PISTONDIR {LEFT,RIGHT};

//SIMPLE WAVE FUNCTIONS

double constant_state_soundspeed(double rho, double pres);

double near_piston_soundspeed(double u1, PISTONDIR dir,
                              double u0, double rho0, double pres0);

double isentropic_relation_density(double a1, double rho0, double pres0);

double isentropic_relation_pressure(double a1, double rho1);

double rarefaction_velocity(double x, double t, PISTONDIR dir,
                            double u0, double a0);

double rarefaction_soundspeed(double x, double t, PISTONDIR dir,
                              double u0, double a0);

//SHOCK WAVE FUNCTIONS

double behind_state_pressure(double rhoa, double pa, double rhob);

double behind_state_specific_volume(double rhoa, double pa, double pb);



#endif
