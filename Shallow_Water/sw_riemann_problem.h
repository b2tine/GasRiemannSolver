#ifndef RIEMANN_PROBLEM_H
#define RIEMANN_PROBLEM_H

#include <secant_method.h>

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

enum class WAVETYPE {SIMPLE,SHOCK,NOWAVE};
enum class DIRECTION {LEFT,RIGHT};


const double G = 9.8; //gravity

struct STATE
{
    double u {0.0};       //velocity
    double h {0.0};       //height
    double p {0.0};       //pressure
    double c {0.0};       //celerity
    
    std::string id;

    STATE() = default;
    
    STATE(const STATE&) = default;
    STATE& operator=(const STATE&) = default;

    STATE(STATE&&) = default;
    STATE& operator=(STATE&&) = default;

    STATE(double U, double H);
    STATE(double U, double H, const std::string& ID);

    void computePressure();
    void computeCelerity();

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


double LeftCenteredWave(double H, STATE* sl, STATE* sl_center);

double RightCenteredWave(double H, STATE* sr_center, STATE* sr);


struct RP_Function
{
    STATE *sleft, *sleft_center, *sright_center, *sright;

    RP_Function() = default;

    //Note: This is not exactly a const method; the STATEs
    //      themselves are being modified in LeftCeteredWave() and
    //      RightCenteredWave() through their pointers.
    //      We need the const modifier in order for RP_Function
    //      to work with the secantMethod() template function.
    double operator()(double H) const
    {
        return LeftCenteredWave(H,sleft,sleft_center)
            - RightCenteredWave(H,sright_center,sright);
    }
};

class RiemannProblem
{
    public:
    
        RiemannProblem(STATE* sL, STATE* sR)
            : sl{sL}, sl_c{new STATE(0,0,"LC")},
            sr_c{new STATE(0,0,"RC")}, sr{sR}
        {
            rpfunc.sleft = sl;
            rpfunc.sleft_center = sl_c;
            rpfunc.sright_center = sr_c;
            rpfunc.sright = sr;
        }

        ~RiemannProblem()
        {
            delete sl_c;
            delete sr_c;
        }

        void solve();
        void printStates();
        void printWaves();

        STATE operator()(double ksi);

        STATE operator()(double x, double t)
        {
            return this->operator()(x/t);
        }

    private:

        RP_Function rpfunc;
        STATE *sl, *sl_c, *sr_c, *sr;

        double H_ctr;
        WAVETYPE LCW {WAVETYPE::NOWAVE};
        WAVETYPE RCW {WAVETYPE::NOWAVE};

        double left_shockspeed;
        double left_trailing_fan_slope;
        double left_leading_fan_slope;
        double right_leading_fan_slope;
        double right_trailing_fan_slope;
        double right_shockspeed;
};


#endif
