#include "riemann_problem.h"
#include <util.h>


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
    std::ifstream infile("in-dambreak");

    std::vector<double> init;
    while (!infile.eof())
    {
        double val;
        infile >> val;
        init.push_back(val);
    }

    double ul = init[0];
    double hl = init[1];

    double ur = init[2];
    double hr = init[3];
    //END INPUT


    STATE sl(ul,hl,"L");   
    STATE sr(ur,hr,"R");   

    RiemannProblem RP(&sl,&sr);
    RP.solve();
    
    //OUTPUT
    RP.printStates();

    std::string outdir("out-dambreak/");
    create_directory(outdir);

    std::ofstream ufile(outdir+"velocity.txt");
    std::ofstream hfile(outdir+"height.txt");
    //std::ofstream afile(outdir+"soundspeed.txt");

    for (int i = 0; i < M; ++i)
    {
        double ksi = xl + i*h;
        if (t > 0.0) ksi /= t;

        STATE U = RP(ksi);
        
        ufile << ksi << " " << U.u << "\n";
        hfile << ksi << " " << U.h << "\n";
        //afile << ksi << " " << U.a << "\n";
    }

    ufile.close();
    hfile.close();
        //afile.close();
    //END OUTPUT

    return 0;
}

