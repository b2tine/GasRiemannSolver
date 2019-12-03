#include "riemann_problem.h"
#include "util.h"


int main()
{
    std::ifstream infile("in-gas");

    std::vector<double> init;
    while (!infile.eof())
    {
        double val;
        infile >> val;
        init.push_back(val);
    }

    double tfinal = 0.1;
    int max_tstep = 2000;
    double default_dt = 0.0005;
    double CFL = 0.75;

    int N = 1000;
    double xmin = 0.0;
    double xmid = 0.5;
    double xmax = 1.0;
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
    std::string outdir("out-gas/");

    std::string density_dir = outdir + "density/";
    create_directory(density_dir);

    std::string velocity_dir = outdir + "velocity/";
    create_directory(velocity_dir);
    
    std::string pressure_dir = outdir + "pressure/";
    create_directory(pressure_dir);


    std::ofstream rhofile(density_dir+"density-0.txt");
    std::ofstream ufile(velocity_dir+"velocity-0.txt");
    std::ofstream pfile(pressure_dir+"pressure-0.txt");

    for (int i = 0; i < N; ++i)
    {
        rhofile << X[i] << " " << U[i].rho << "\n";
        ufile << X[i] << " " << U[i].u << "\n";
        pfile << X[i] << " " << U[i].p << "\n";
    }

    rhofile.close();
    ufile.close();
    pfile.close();

    //Start Up Step
    //TODO: are we doing this correctly, and is it neccesary?
    RiemannProblem RP_StartUp(&ULeftDirichlet,&URightDirichlet);
    RP_StartUp.solve();
    STATE V_StartUp = RP_StartUp(0.0);

    double u_start = fabs(V_StartUp.u);
    double a_start = V_StartUp.a;

    double max_speed =
        std::max((std::max(u_start,fabs(u_start-a_start))),(fabs(u_start+a_start)));
    double max_dt = CFL*dx/max_speed;

    double time = 0.0;
    std::ofstream logfile(outdir + "log.txt");

    //Time Marching
    for (int ts = 1; ts <= max_tstep; ++ts)
    {
        max_speed = 0.0;

        double dt = std::min(default_dt,max_dt);
        time += dt;

        //Domain Interior
        for (int i = 1; i < N-1; ++i)
        {
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

            //record max_speed for new dt
            double u = fabs(Unew[i].u);
            double a = Unew[i].a;
            double curr_max_speed = std::max((std::max(u,fabs(u-a))),(fabs(u+a)));
            if (max_speed < curr_max_speed)
                max_speed = curr_max_speed;
        }

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

        //write output files
        std::string tstep = std::to_string(ts);
        rhofile.open(density_dir+"density-"+tstep+".txt");
        ufile.open(velocity_dir+"velocity-"+tstep+".txt");
        pfile.open(pressure_dir+"pressure-"+tstep+".txt");

        for (int i = 0; i < N; ++i)
        {
            rhofile << X[i] << " " << U[i].rho << "\n";
            ufile << X[i] << " " << U[i].u << "\n";
            pfile << X[i] << " " << U[i].p << "\n";
        }

        rhofile.close();
        ufile.close();
        pfile.close();

        logfile << "step = " << ts << " ";
        logfile << "dt = " << dt << " ";
        logfile << "time = " << time << "\n\n";

        if (time >= tfinal)
            break;
    }

    logfile.close();

    return 0;
}

