#include "sw_riemann_problem.h"
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

    //TODO: read these from input file
    double tfinal = 0.01;
    int max_tstep = 2000;
    double default_dt = 0.0005;
    double CFL = 0.75;

    int N = 1000;
    double xmin = 0.0;
    double xmid = 0.5;
    double xmax = 1.0;
    double dx = (xmax - xmin)/((double)N);

    double X[N];

    double Q[N][2];
    double Qnew[N][2];
    
    STATE U[N];
    STATE Unew[N];

    //generate grid and set initial condition
    for (int i = 0; i < N; ++i)
    {
        double velo, height;

        X[i] = xmin + ((double) i + 0.5)*dx;

        if (X[i] < xmid)
        {
            velo = init[0];
            height = init[1];
        }
        else
        {
            velo = init[2];
            height = init[3];
        }

        Q[i][0] = height;
        Q[i][1] = velo*height;

        U[i].u = velo;
        U[i].h = height;
    }
    
    //Dirichlet Boundaries
    double QLeftDirichlet[2] = {Q[0][0], Q[0][1]};
    STATE ULeftDirichlet = U[0];

    double QRightDirichlet[2] = {Q[N-1][0], Q[N-1][1]};
    STATE URightDirichlet = U[N-1];


    //write initial condition output files
    std::string outdir(out_name + "/");

    std::string velocity_dir = outdir + "velocity/";
    create_directory(velocity_dir);
    
    std::string height_dir = outdir + "height/";
    create_directory(height_dir);

    std::string pressure_dir = outdir + "pressure/";
    create_directory(pressure_dir);

    std::string celerity_dir = outdir + "celerity/";
    create_directory(celerity_dir);

    char ufilename[100];
    sprintf(ufilename,"%s/velocity-%04d.txt",velocity_dir.c_str(),0);
    std::ofstream ufile(ufilename);

    char hfilename[100];
    sprintf(hfilename,"%s/height-%04d.txt",height_dir.c_str(),0);
    std::ofstream hfile(hfilename);

    char pfilename[100];
    sprintf(pfilename,"%s/pressure-%04d.txt",pressure_dir.c_str(),0);
    std::ofstream pfile(pfilename);

    char cfilename[100];
    sprintf(cfilename,"%s/celerity-%04d.txt",celerity_dir.c_str(),0);
    std::ofstream cfile(cfilename);
    
    for (int i = 0; i < N; ++i)
    {
        ufile << X[i] << " " << U[i].u << "\n";
        hfile << X[i] << " " << U[i].h << "\n";
        pfile << X[i] << " " << U[i].p << "\n";
        cfile << X[i] << " " << U[i].c << "\n";
    }

    ufile.close();
    hfile.close();
    pfile.close();
    cfile.close();

    //Start Up Step
    RiemannProblem RP_StartUp(&ULeftDirichlet,&URightDirichlet);
    RP_StartUp.solve();
    STATE V_StartUp = RP_StartUp(0.0);

    double u_start = fabs(V_StartUp.u);
    double c_start = V_StartUp.c;

    double max_speed =
        std::max((std::max(u_start,fabs(u_start-c_start))),
                 (fabs(u_start+c_start)));
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
            QFlux[0] = VP.u*VP.h - VM.u*VM.h;
            QFlux[1] = VP.u*VP.u*VP.h + VP.p - (VM.u*VM.u*VM.h + VM.p);

            //Update Conserved Variables and Field States
            double height = Q[i][0] - QFlux[0]*dt/dx;
            double velheight = Q[i][1] - QFlux[1]*dt/dx;

            Qnew[i][0] = height;
            Qnew[i][1] = velheight;

            Unew[i].h = height;
            Unew[i].u = velheight/height;
            Unew[i].computePressure();
            Unew[i].computeCelerity();

            //Record max_speed next time step
            double u = fabs(Unew[i].u);
            double c = Unew[i].c;
            double curr_max_speed = std::max((std::max(u,fabs(u-c))),(fabs(u+c)));
            if (max_speed < curr_max_speed)
                max_speed = curr_max_speed;
        }

        //Update time step data
        if (max_speed > 0.0)
            max_dt = CFL*dx/max_speed;
 
        //Dirichlet Boundary Conditions
        Qnew[0][0] = QLeftDirichlet[0];
        Qnew[0][1] = QLeftDirichlet[1];
        Unew[0] = ULeftDirichlet;

        Qnew[N-1][0] = QRightDirichlet[0];
        Qnew[N-1][1] = QRightDirichlet[1];
        Unew[N-1] = URightDirichlet;

        //copy new into old
        for (int i = 0; i < N; ++i)
        {
            Q[i][0] = Qnew[i][0];
            Q[i][1] = Qnew[i][1];
            U[i] = Unew[i];
        }

        //write output files
        sprintf(ufilename,"%s/velocity-%04d.txt",velocity_dir.c_str(),ts);
        sprintf(hfilename,"%s/height-%04d.txt",height_dir.c_str(),ts);
        sprintf(pfilename,"%s/pressure-%04d.txt",pressure_dir.c_str(),ts);
        sprintf(cfilename,"%s/celerity-%04d.txt",celerity_dir.c_str(),ts);
        
        ufile.open(ufilename);
        hfile.open(hfilename);
        pfile.open(pfilename);
        cfile.open(cfilename);

        for (int i = 0; i < N; ++i)
        {
            ufile << X[i] << " " << U[i].u << "\n";
            hfile << X[i] << " " << U[i].h << "\n";
            pfile << X[i] << " " << U[i].p << "\n";
            cfile << X[i] << " " << U[i].c << "\n";
        }

        ufile.close();
        hfile.close();
        pfile.close();
        cfile.close();

        logfile << "step = " << ts << "   ";
        logfile << "time = " << time << "   ";
        logfile << "dt = " << dt << "\n\n";

        if (time >= tfinal)
            break;
    }

    logfile.close();

    return 0;
}

