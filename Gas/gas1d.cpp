#include "riemann_problem.h"
#include <util.h>


int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        printf("ERROR: Require the input file name \
                and output directory name.\n");
        exit(1);
    }

    std::string in_name = argv[1];
    std::string out_name = argv[2];

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

    double tfinal = init[6];
    int max_tstep = 10000;
    double default_dt = 0.0005;
    double CFL = 0.75;

    int N = 1000;
    double xmin = init[7];
    double xmid = init[8];
    double xmax = init[9];
    double dx = (xmax - xmin)/((double)N);

    double X[N];

    double Q[N][3];
    double Qnew[N][3];
    
    STATE U[N];
    STATE Unew[N];

    //generate grid and set initial condition
    for (int i = 0; i < N; ++i)
    {
        double dens, velo, pres;

        X[i] = xmin + ((double) i + 0.5)*dx;

        if (X[i] < xmid)
        {
            dens = init[0];
            velo = init[1];
            pres = init[2];
        }
        else
        {
            dens = init[3];
            velo = init[4];
            pres = init[5];
        }

        Q[i][0] = dens;
        Q[i][1] = dens*velo;
        Q[i][2] = 0.5*dens*velo*velo + pres/(GAMMA - 1.0);

        U[i].rho = dens;
        U[i].u = velo;
        U[i].p = pres;
        U[i].computeSoundSpeed();
    }
    
    //Dirichlet Boundaries
    double QLeftDirichlet[3] = {Q[0][0], Q[0][1], Q[0][2]};
    STATE ULeftDirichlet = U[0];

    double QRightDirichlet[3] = {Q[N-1][0], Q[N-1][1], Q[N-1][2]};
    STATE URightDirichlet = U[N-1];


    //write initial condition output files
    std::string outdir(out_name + "/");

    std::string density_dir = outdir + "density/";
    create_directory(density_dir);

    std::string velocity_dir = outdir + "velocity/";
    create_directory(velocity_dir);
    
    std::string pressure_dir = outdir + "pressure/";
    create_directory(pressure_dir);

    std::string soundspeed_dir = outdir + "soundspeed/";
    create_directory(soundspeed_dir);

    char rhofilename[100];
    sprintf(rhofilename,"%s/density-%04d.txt",density_dir.c_str(),0);
    std::ofstream rhofile(rhofilename);

    char ufilename[100];
    sprintf(ufilename,"%s/velocity-%04d.txt",velocity_dir.c_str(),0);
    std::ofstream ufile(ufilename);

    char pfilename[100];
    sprintf(pfilename,"%s/pressure-%04d.txt",pressure_dir.c_str(),0);
    std::ofstream pfile(pfilename);

    char afilename[100];
    sprintf(afilename,"%s/soundspeed-%04d.txt",soundspeed_dir.c_str(),0);
    std::ofstream afile(afilename);
    
    for (int i = 0; i < N; ++i)
    {
        rhofile << X[i] << " " << U[i].rho << "\n";
        ufile << X[i] << " " << U[i].u << "\n";
        pfile << X[i] << " " << U[i].p << "\n";
        afile << X[i] << " " << U[i].a << "\n";
    }

    rhofile.close();
    ufile.close();
    pfile.close();
    afile.close();

    //Start Up Step
    RiemannProblem RP_StartUp(&ULeftDirichlet,&URightDirichlet);
    RP_StartUp.solve();
    STATE V_StartUp = RP_StartUp(0.0);

    double u_start = fabs(V_StartUp.u);
    double a_start = V_StartUp.a;

    double max_speed =
        std::max((std::max(u_start,fabs(u_start-a_start))),(fabs(u_start+a_start)));
    double max_dt = CFL*dx/max_speed;

    //Time Marching
    double time = 0.0;
    std::ofstream logfile(outdir + "log.txt");

    for (int ts = 1; ts <= max_tstep; ++ts)
    {
        max_speed = 0.0;

        double dt = std::min(default_dt,max_dt);
        time += dt;

        //Domain Interior
        for (int i = 1; i < N-1; ++i)
        {
            //Compute Godunov Flux
            RiemannProblem RP_plus(&U[i],&U[i+1]);
            RP_plus.solve();
            STATE VP = RP_plus(0.0);

            RiemannProblem RP_minus(&U[i-1],&U[i]);
            RP_minus.solve();
            STATE VM = RP_minus(0.0);

            double QFlux[3];
            QFlux[0] = VP.rho*VP.u - VM.rho*VM.u;
            QFlux[1] = VP.rho*VP.u*VP.u + VP.p - (VM.rho*VM.u*VM.u + VM.p);
            QFlux[2] = VP.u*(0.5*VP.rho*VP.u*VP.u + VP.p/(GAMMA - 1.0) + VP.p)
                        - VM.u*(0.5*VM.rho*VM.u*VM.u + VM.p/(GAMMA - 1.0) + VM.p);

            //Update Conserved Variables and Field States
            double dens = Q[i][0] - QFlux[0]*dt/dx;
            double momn = Q[i][1] - QFlux[1]*dt/dx;
            double energy = Q[i][2] - QFlux[2]*dt/dx;

            Qnew[i][0] = dens;
            Qnew[i][1] = momn;
            Qnew[i][2] = energy;

            Unew[i].rho = dens;
            Unew[i].u = momn/dens;
            Unew[i].p = (energy - 0.5*momn*momn/dens)*(GAMMA - 1.0);
            Unew[i].computeSoundSpeed();

            //Record max_speed next time step
            double u = fabs(Unew[i].u);
            double a = Unew[i].a;
            double curr_max_speed = std::max((std::max(u,fabs(u-a))),(fabs(u+a)));
            if (max_speed < curr_max_speed)
                max_speed = curr_max_speed;
        }

        //Update time step data
        if (max_speed > 0.0)
            max_dt = CFL*dx/max_speed;
 
        //Dirichlet Boundary Conditions
        Qnew[0][0] = QLeftDirichlet[0];
        Qnew[0][1] = QLeftDirichlet[1];
        Qnew[0][2] = QLeftDirichlet[2];
        Unew[0] = ULeftDirichlet;

        Qnew[N-1][0] = QRightDirichlet[0];
        Qnew[N-1][1] = QRightDirichlet[1];
        Qnew[N-1][2] = QRightDirichlet[2];
        Unew[N-1] = URightDirichlet;

        //copy new into old
        for (int i = 0; i < N; ++i)
        {
            Q[i][0] = Qnew[i][0];
            Q[i][1] = Qnew[i][1];
            Q[i][2] = Qnew[i][2];
            U[i] = Unew[i];
        }

        if (ts % 5 == 0)
        {
            //write output files
            sprintf(rhofilename,"%s/density-%04d.txt",density_dir.c_str(),ts);
            sprintf(ufilename,"%s/velocity-%04d.txt",velocity_dir.c_str(),ts);
            sprintf(pfilename,"%s/pressure-%04d.txt",pressure_dir.c_str(),ts);
            sprintf(afilename,"%s/soundspeed-%04d.txt",soundspeed_dir.c_str(),ts);
            
            rhofile.open(rhofilename);
            ufile.open(ufilename);
            pfile.open(pfilename);
            afile.open(afilename);

            for (int i = 0; i < N; ++i)
            {
                rhofile << X[i] << " " << U[i].rho << "\n";
                ufile << X[i] << " " << U[i].u << "\n";
                pfile << X[i] << " " << U[i].p << "\n";
                afile << X[i] << " " << U[i].a << "\n";
            }

            rhofile.close();
            ufile.close();
            pfile.close();
            afile.close();
        }

        logfile << "step = " << ts << "   ";
        logfile << "time = " << time << "   ";
        logfile << "dt = " << dt << "\n\n";

        if (time >= tfinal)
            break;
    }

    logfile.close();

    return 0;
}

