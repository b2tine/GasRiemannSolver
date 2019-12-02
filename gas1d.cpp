#include "riemann_problem.h"
#include "util.h"


int main()
{
    //INPUT: Left and Right states, sl and sr

    std::ifstream infile("in-gas");

    std::vector<double> init;
    while (!infile.eof())
    {
        double val;
        infile >> val;
        init.push_back(val);
    }

    //1. generate initial condition and grid:
    //   4 arrays, domain, vel, density, pressure

    int N = 100;
    int max_tstep = 10;

    double xmin = 0.0;
    double xmid = 0.5;
    double xmax = 1.0;
    double dx = (xmax - xmin)/((double)N-1);

    double X[N];
    STATE U[N];

    for (int i = 0; i < N; ++i)
    {
        double dens, vel, pres;

        X[i] = xmin + ((double) i + 0.5)*dx;

        if (X[i] < xmid)
        {
            dens = init[0];
            vel = init[1];
            pres = init[2];
        }
        else
        {
            dens = init[3];
            vel = init[4];
            pres = init[5];
        }

        U[i].rho = dens;
        U[i].u = vel;
        U[i].p = pres;
    }


    for (int ts = 1; ts <= max_tstep; ++ts)
    {
        //interior
        for (int i = 1; i < N-1; ++i)
        {
            //TODO: these center states should be made data
            //      members of the RiemannProblem class.
            STATE slc_left(0,0,0);
            STATE src_left(0,0,0);

            RiemannProblem RP_left(&U[i-1],&slc_left,&src_left,&U[i]);
            RP_left.solve();

            STATE slc_right(0,0,0);
            STATE src_right(0,0,0);

            RiemannProblem RP_right(&U[i],&slc_right,&src_right,&U[i+1]);
            RP_right.solve();

            //un[i] = u[i] - ( GFlux(f,u[i],u[i+1]) - GFlux(f,u[i-1],u[i]) )*dt/dx
    
        }

        //periodic boundary
        
            //un[0] = un[N-2];
            //un[N-1] = un[1];
    
    }

    return 0;
}
