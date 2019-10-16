#include "riemann_problem.h"
#include "util.h"


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

    
    //INPUT: Left and Right states, sl and sr
    std::ifstream infile("in-RP");

    std::vector<double> init;
    while (!infile.eof())
    {
        double val;
        infile >> val;
        init.push_back(val);
    }

    double ul = init[0];
    double rhol = init[1];
    double pl = init[2];

    double ur = init[3];
    double rhor = init[4];
    double pr = init[5];
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
    std::string outdir("out-RP/");
    create_directory(outdir);

    std::ofstream ufile(outdir+"velocity.txt");
    std::ofstream rhofile(outdir+"density.txt");
    std::ofstream pfile(outdir+"pressure.txt");
    std::ofstream afile(outdir+"soundspeed.txt");

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

