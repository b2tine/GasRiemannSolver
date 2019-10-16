#include "riemann_problem.h"


int main(int argc, char* argv[])
{
    double xl = atof(argv[1]);
    double xr = atof(argv[2]);
    double t = atof(argv[3]);
    //double m = atof(argv[3]);
    //double t = atof(argv[4]);

    assert (xl <= xr);
    assert (t > 0.0);

    //GIVEN: left and right states, sl and sr
    double ul = 0.0;
    double rhol = 1.225;
    double pl = 100000;
    double al = constant_state_soundspeed(rhol,pl);

    double ur = -1.0;
    double rhor = 1.1;
    double pr = 85000;
    double ar = constant_state_soundspeed(rhor,pr);

    STATE sl = {ul, rhol, pl, al, "L"};      //sl.id = "L";
    STATE sl_c = {0, 0, 0, 0, "LC"};         //sl_c.id = "LC"
    STATE sr_c = {0, 0, 0, 0, "RC"};         //sr_c.id = "RC"
    STATE sr = {ur, rhor, pr, ar, "R"};      //sr.id = "R";

    RiemannProblem RP(&sl,&sl_c,&sr_c,&sr);
    
    std::cout << sl << "\n";
    std::cout << sl_c << "\n";
    std::cout << sr_c << "\n";
    std::cout << sr << "\n\n";

    /*
    //sample (x,t) point 
    double x = 0.5;
    double t = 0.02;

    double ksi = x/t;
    */

    int M = 25;
    double h = (xr-xl)/M;
    for (int i = 0; i < M+1; ++i)
    {
        double x = xl + i*h;
        double ksi = x/t;
        double u_riemann = RP(ksi);
        printf("u_riemann(%g) = %g\n",ksi,u_riemann);
    }

    //printf("u_riemann(%g/%g = %g) = %g\n",x,t,ksi,u_riemann);

    return 0;
}

