#include "sw_riemann_problem.h"
#include <util.h>


int main(int argc, char* argv[])
{
    //TODO: read these from input file
    /*
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
    */

    //INPUT: Left and Right states, sl and sr
    if (argc < 3)
    {
        printf("ERROR: Require the input file name \
                and output directory name.\n");
        exit(1);
    }

    std::string in_name = argv[1];
    std::string out_name = argv[2];

    //Read input file
    std::ifstream infile(in_name);

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

    //Initalize states and solve Riemann Problem
    STATE sl(ul,hl,"L");   
    STATE sr(ur,hr,"R");   

    RiemannProblem RP(&sl,&sr);
    RP.solve();
    RP.printStates();
    RP.printWaves();
    
    //Write output files
    std::string outdir(out_name + "/");
    create_directory(outdir);

    std::ofstream ufile(outdir+"velocity.txt");
    std::ofstream hfile(outdir+"height.txt");
    //std::ofstream afile(outdir+"soundspeed.txt");

    int M = 500;
    double xl = -15;
    double xr = 15;
    double h = (xr-xl)/M;
    double t = 1.0;

    for (int i = 0; i < M; ++i)
    {
        double ksi = xl + i*h;
        ksi /= t;

        STATE U = RP(ksi);
        
        ufile << ksi << " " << U.u << "\n";
        hfile << ksi << " " << U.h << "\n";
        //afile << ksi << " " << U.a << "\n";
    }

    ufile.close();
    hfile.close();
        //afile.close();

    return 0;
}

