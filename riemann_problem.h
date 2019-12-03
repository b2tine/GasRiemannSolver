#ifndef RIEMANN_PROBLEM_H
#define RIEMANN_PROBLEM_H

#include "secant_method.h"

#include <iostream>
#include <exception>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <limits>


const double HUGE = std::numeric_limits<double>::max();


const double GAMMA = 1.4;

enum class WAVETYPE {SHOCK,SIMPLE};
enum class DIRECTION {LEFT,RIGHT};

struct STATE
{
    double rho {-1.0};    //density
    double u {0.0};       //velocity
    double p {-1.0};      //pressure
    double a {-1.0};      //soundspeed
    //double gamma;       //specific heat ratio
    
    std::string id;

    STATE() = default;
    
    STATE(const STATE&) = default;
    STATE& operator=(const STATE&) = default;

    STATE(STATE&&) = default;
    STATE& operator=(STATE&&) = default;

    STATE(double RHO, double U, double P);
    STATE(double RHO, double U, double P, const std::string& ID);
    STATE(double RHO, double U, double P, double A);
    STATE(double RHO, double U, double P, double A, const std::string& ID);

    std::string printinfo() const;

    void print() const
    {
        printf("%s",printinfo().c_str());
    }

    friend std::ostream& operator << (std::ostream& os, const STATE& s)
    {
        os << s.printinfo();
        return os;
    }
};


double LeftCenteredWave(double Pslip, STATE* sl, STATE* sl_center);

double RightCenteredWave(double Pslip, STATE* sr_center, STATE* sr);


struct RP_Function
{
    STATE *sleft, *sleft_center, *sright_center, *sright;

    RP_Function() = default;

    RP_Function(STATE* sL, STATE* sLC, STATE* sRC, STATE* sR)
        : sleft{sL}, sleft_center{sLC}, sright_center{sRC}, sright{sR}
    {}

    //Note: This is not exactly a const method; the STATEs
    //      themselves are being modified in LeftCeteredWave() and
    //      RightCenteredWave() through their pointers.
    //      We need the const modifier in order for RP_Function
    //      to work with the secantMethod() template function.
    double operator()(double P) const
    {
        return LeftCenteredWave(P,sleft,sleft_center)
            - RightCenteredWave(P,sright_center,sright);
    }
};

class RiemannProblem
{
    public:
    
        RiemannProblem(STATE* sL, STATE* sR)
            : sl{sL}, sl_c{new STATE(0,0,0,"LC")},
            sr_c{new STATE(0,0,0,"RC")}, sr{sR}
        {
            rpfunc.sleft = sL;
            rpfunc.sleft_center = sl_c;
            rpfunc.sright_center = sr_c;
            rpfunc.sright = sR;
        }

        ~RiemannProblem()
        {
            delete sl_c;
            delete sr_c;
        }

        void solve();

        STATE operator()(double ksi);

        STATE operator()(double x, double t)
        {
            return this->operator()(x/t);
        }

        //TODO: printing functions

        void printStates();


    private:

        STATE *sl, *sl_c, *sr_c, *sr;
        RP_Function rpfunc;

        double Pslip;
        WAVETYPE LCW, RCW;

        double left_shockspeed {HUGE};
        double left_trailing_fan_slope {HUGE};
        double left_leading_fan_slope {HUGE};
        double slip_slope {0.0};
        double right_trailing_fan_slope {-HUGE};
        double right_leading_fan_slope {-HUGE};
        double right_shockspeed {-HUGE};

        void detectVacuumState();
};


class VacuumStateException : public std::runtime_error
{
    public:

        VacuumStateException(const std::string& arg)
            : std::runtime_error{arg}
        {
            std::ostringstream oss;
            oss << "ERROR: Vacuum State detected " << arg << ".\n";
            message = oss.str();
        }

        VacuumStateException(const std::string& arg,
                             const std::vector<STATE*>& vstates)
            : VacuumStateException{arg}
        {
            std::ostringstream oss;
            for (STATE* s : vstates)
                oss << *s << "\n";
            message += oss.str();
        }

        const char* what() const noexcept override
        {
            return message.c_str();
        }

    private:

        std::string message;
};


//SHOCK WAVE FUNCTIONS

double behind_state_specific_volume(double rhoa, double pa, double pb);

double behind_state_pressure(double rhoa, double pa, double rhob);


//SIMPLE WAVE FUNCTIONS

double constant_state_soundspeed(double rho, double pres);

double near_piston_soundspeed(double u1, DIRECTION dir,
                              double u0, double rho0, double pres0);

double isentropic_relation_density(double a1, double rho0, double pres0);

double isentropic_relation_pressure(double a1, double rho1);

double rarefaction_velocity(double ksi, DIRECTION dir, double u0, double a0);

double rarefaction_velocity_xt(double x, double t,
                            DIRECTION dir, double u0, double a0);

double rarefaction_soundspeed(double ksi, DIRECTION dir, double u0, double a0);

double rarefaction_soundspeed_xt(double x, double t,
                              DIRECTION dir, double u0, double a0);


#endif
