#include "riemann_problem.h"

double pressure_behind(double rhoa, double pa, double rhob);
double specific_volume_behind(double rhoa, double pa, double pb);

//LFS: ul, ur > 0
//     u0 = ul = ua and u1 = ur = ub

int main(int argc, char* argv[])
{
    //double gamma = 1.4;

    
    //0.given 4 inputs: ul, rhol, pl and either rhor OR pr
    //0.given 4 inputs: u0, rho0, p0 and either rho1 OR p1
    //0.given 4 inputs: ua, rhoa, pa and either rhob OR pb
    
    double ul = 1.0;            //ua
    double rhol = 1.0;          //rhoa
    double taul = 1.0/rhol;     //taua
    double pl = 900000;         //pa

    //TODO: test giving rhob instead of pb (rhor instead of pr)
    //double rhor = 1.225;      //rhob
    //double taur = 1.0/rhor;   //taub
    double pr = 100000;         //pb

    STATE sl = {ul, rhol, pl};  //ahead state, sa
    //STATE sr = {ur, rhor, pr};

    //1.use Hugoniot function to compute pb OR rhob (whichever not given)
    //   --since pr given, compute rhor (or taur if easier)

    //double pr = pressure_behind(rhol,pl,rhor);        //pb
    double taur = specific_volume_behind(rhol,pl,pr);   //taub
    double rhor = 1.0/taur;                             //rhob
    
    //NOTE: M is mass flux in the shock-stationary frame, not the mach number

    //2.compute M from M^2 = (pr - pl)/(taul - taur)
    //2.compute M from M^2 = (p1 - p0)/(tau0 - tau1)
    //2.compute M from M^2 = (pb - pa)/(taua - taub)

    double M = std::sqrt((pr - pl)/(taul - taur));

    //3.set M = (pr - pl)/(ul - ur) and compute ur
    //3.set M = (p1 - p0)/(u0 - u1) and compute u1
    //3.set M = (pb - pa)/(ua - ub) and compute ub

    double ur = ul - (pr - pl)/M;

    //4. check that entropy condition satisfied;
    //  i.e. rhoa < rhob && pa < pb

    bool entropy_condition = true;
    if (rhol >= rhor)
    {
        entropy_condition = false;
        printf("ahead state density not less than behind state density\n");
    }
    if (pl >= pr)
    {
        entropy_condition = false;
        printf("ahead state pressure not less than behind state pressure\n");
    }
    if (!entropy_condition)
    {
        printf("entropy condition not satisfied\n");
    }

    return 0;
}

double pressure_behind(double rhoa, double pa, double rhob)
{
    double GP = GAMMA + 1.0;
    double GM = GAMMA - 1.0;
    return (GP/rhoa - GM/rhob)/(GP/rhob - GM/rhoa)*pa;
}

double specific_volume_behind(rhoa,pa,pb)
{
    double GP = GAMMA + 1.0;
    double GM = GAMMA - 1.0;
    return (GP*pa + GM*pb)/(GP*pb + GM*pa)/rhoa;
}

