#include "riemann_problem.h"
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
    
    std::string in_name(argv[1]);
    std::string out_name(argv[2]);

    //Read input file
    std::vector<double> init;
    std::ifstream infile(in_name);
    while (!infile.eof())
    {
        double val;
        infile >> val;
        init.push_back(val);
    }
    infile.close();

    double rhol = init[0];
    double ul = init[1];
    double pl = init[2];

    double rhor = init[3];
    double ur = init[4];
    double pr = init[5];

    //Initalize states and solve Riemann Problem
    STATE sl(rhol,ul,pl,"L");   
    STATE sr(rhor,ur,pr,"R");   

    RiemannProblem RP(&sl,&sr);
    RP.solve();
    
    RP.printStates();
    RP.printWaves();
    
    //Write output files
    std::string outdir(out_name + "/");
    create_directory(outdir);

    std::ofstream rhofile(outdir+"density.txt");
    std::ofstream ufile(outdir+"velocity.txt");
    std::ofstream pfile(outdir+"pressure.txt");
    std::ofstream afile(outdir+"soundspeed.txt");

    int M = 400;
    double xl = -1.0;
    double xr = 1.0;
    double h = (xr-xl)/M;
    double t = 0.4;

    for (int i = 0; i < M; ++i)
    {
        double ksi = xl + i*h;
        ksi /= t;

        STATE uR = RP(ksi);
        
        rhofile << ksi << " " << uR.rho << "\n";
        ufile << ksi << " " << uR.u << "\n";
        pfile << ksi << " " << uR.p << "\n";
        afile << ksi << " " << uR.a << "\n";
    }

    rhofile.close();
    ufile.close();
    pfile.close();
    afile.close();

    return 0;
}

