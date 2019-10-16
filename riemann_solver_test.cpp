#include "riemann_problem.h"


int main(int argc, char* argv[])
{
    //INIT OUTPUT:
    double xl = atof(argv[1]);
    double xr = atof(argv[2]);
    int M = atoi(argv[3]);
    assert (xl <= xr && M >= 1);
    
    double t = -1.0;
    if (argc > 4)
    {
        t = atof(argv[4]);
        assert (t > 0.0);
    }

    double h = (M == 1) ? xr-xl : (xr-xl)/(M-1);
    if (h == 0.0) M = 1;
    //END INIT OUTPUT

    
    //TODO: read from input file
    //INPUT: left and right states, sl and sr
    double ul = 0.0;
    double rhol = 1.225;
    double pl = 100000;

    double ur = -1.0;
    double rhor = 1.1;
    double pr = 85000;
    //END INPUT


    STATE sl(ul,rhol,pl,"L");   
    STATE sl_c(0,0,0,"LC");   
    STATE sr_c(0,0,0,"RC");   
    STATE sr(ur,rhor,pr,"R");   

    RiemannProblem RP(&sl,&sl_c,&sr_c,&sr);
    
    std::cout << sl << "\n";
    std::cout << sl_c << "\n";
    std::cout << sr_c << "\n";
    std::cout << sr << "\n";

    
    //OUTPUT
    std::ofstream ufile("velocity-xt.txt");
    std::ofstream rhofile("density-xt.txt");
    std::ofstream pfile("pressure-xt.txt");
    std::ofstream afile("soundspeed-xt.txt");

    for (int i = 0; i < M; ++i)
    {
        double ksi = xl + i*h;
        if (t > 0.0) ksi /= t;

        STATE uR = RP(ksi);
        
        ufile << ksi << " " << uR.u << "\n";
        rhofile << ksi << " " << uR.rho << "\n";
        pfile << ksi << " " << uR.p << "\n";
        afile << ksi << " " << uR.a << "\n";
    }

    ufile.close();
    rhofile.close();
    pfile.close();
    afile.close();
    //END OUTPUT

    return 0;
}

