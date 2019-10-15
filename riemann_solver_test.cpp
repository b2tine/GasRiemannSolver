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
    //TODO timings for othe guess (comment block below) instead
    //of just using pl and pr.
    
    /* 
        double frac_prange = 0.05*std::abs(pl-pr);
        double pg0 = 0.5*(pl+pr);
        double pg1 = pg0 + randsign()*frac_prange;
        printf("frac_prange = %g\n",frac_prange);
        printf("pg0 = %g\npg1 = %g\n\n",pg0,pg1);
    */

    //Solve F(Pslip) = ul_star(P_slip) - ur-star(P_slip) = 0
    //using the secant method to find the pressure at the contact discontinuity
    
    RiemannProblem RP(&sl,&sl_c,&sr_c,&sr);
    //RP_Function F(&sl,&sl_c,&sr_c,&sr);
    //double Pslip = secantMethod(F,pl,pr);
        //double Pslip = secantMethod(F,pg0,pg1);
    
    //std::cout << "Pslip = " << Pslip << "\n";

    std::cout << "sl = " << sl;
    std::cout << "sl_c = " << sl_c;
    std::cout << "sr_c = " << sr_c;
    std::cout << "sr = " << sr;

    
    //sample (x,t) point 
    double x = 0.5;
    double t = 0.02;

    double ksi = x/t;
    double u_riemann = RP(ksi);

    std::cout << "u_riemann = " << u_riemann << "\n";

    return 0;
}

