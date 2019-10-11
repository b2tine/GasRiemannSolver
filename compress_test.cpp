#include "riemann_problem.h"

//const double GAMMA = 1.4;

double behind_state_pressure(double rhoa, double pa, double rhob);
double behind_state_specific_volume(double rhoa, double pa, double pb);

//LFS: ul, ur > 0
//     u0 = ul = ua and u1 = ur = ub

int main(int argc, char* argv[])
{
    //double gamma = 1.4;

    //TODO: generalize to ahead and behind state variables
    
    //0.given 4 inputs: ul, rhol, pl and either rhor OR pr
    //0.given 4 inputs: ua, rhoa, pa and either rhob OR pb
    
    double ul = 360.0;            //ua
    double rhol = 1.0;          //rhoa
    double taul = 1.0/rhol;     //taua
    double pl = 90000;         //pa
    double al = constant_state_soundspeed(rhol,pl);

    //NOTE: ul must be greater than al for RIGHT piston compression.
    //      more generally, ua most be supersonic.

    printf("ul = %g\n",ul);
    printf("rhol = %g\n",rhol);
    printf("pl = %g\n",pl);
    printf("al = %g\n\n",al);

    //TODO: check STATE struct to see which variable has been provided?
    double rhor = 1.225;      //rhob
    double taur = 1.0/rhor;   //taub
    //double pr = 100000;         //pb

    //STATE sl = {ul, rhol, pl};  //ahead state, sa
    //STATE sr = {ur, rhor, pr};

    //1.use Hugoniot function to compute pb OR rhob (whichever not given)
    //   --since pr given, compute rhor (or taur if easier)

    double pr = behind_state_pressure(rhol,pl,rhor);            //pb
    //double taur = behind_state_specific_volume(rhol,pl,pr);   //taub
    //double rhor = 1.0/taur;                                   //rhob
    double ar = constant_state_soundspeed(rhor,pr);

    //NOTE: M is mass flux in the shock-stationary frame, not the mach number

    //2.compute M from M^2 = (pr - pl)/(taul - taur)
    //2.compute M from M^2 = (pb - pa)/(taua - taub)

    double M = std::sqrt((pr - pl)/(taul - taur));

    //3.set M = (pr - pl)/(ul - ur) and compute ur
    //3.set M = (pb - pa)/(ua - ub) and compute ub

    double ur = ul - (pr - pl)/M;

    printf("ur = %g\n",ur);
    printf("rhor = %g\n",rhor);
    printf("pr = %g\n",pr);
    printf("ar = %g\n\n",ar);

    //TODO: check prandtl relation, ua and ub have same sign for shock

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

    double shock_speed = (rhol*ul - rhor*ur)/(rhol - rhor);
    printf("shock_speed = %g\n",shock_speed);

    return 0;
}

double behind_state_pressure(double rhoa, double pa, double rhob)
{
    double GP = GAMMA + 1.0;
    double GM = GAMMA - 1.0;
    return (GP/rhoa - GM/rhob)/(GP/rhob - GM/rhoa)*pa;
}

double behind_state_specific_volume(double rhoa, double pa, double pb)
{
    double GP = GAMMA + 1.0;
    double GM = GAMMA - 1.0;
    return (GP*pa + GM*pb)/(GP*pb + GM*pa)/rhoa;
}

