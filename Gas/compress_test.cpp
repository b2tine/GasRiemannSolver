#include "riemann_problem.h"

//const double GAMMA = 1.4;

//RIGHT piston compression
//LFS: ul, ur > 0
//     ul = ua and ur = ub

//LEFT piston compression
//RFS: ul, ur < 0
//     ul = ub and ur = ua

int main(int argc, char* argv[])
{
    double sign = 1.0;
    DIRECTION dir = DIRECTION::RIGHT;
    if (argc > 1)
    {
        if (argv[1][0] == 'l' || argv[1][0] == 'L')
        {
            sign = -1.0;
            dir = DIRECTION::LEFT;
        }
    }


    //TODO: generalize to ahead and behind state variables
    
    //0.given 4 inputs: ul, rhol, pl and either rhor OR pr
    //0.given 4 inputs: ua, rhoa, pa and either rhob OR pb
    
    double Uleft = 385.0;        // m/s
    double ua = sign*Uleft;

    printf("%s piston compression speed U = %g m/s\n\n",
            (dir == DIRECTION::LEFT) ? "LEFT" : "RIGHT", ua);

    double rhoa = 1.0;          // kg/m^3
    double taua = 1.0/rhoa;     // m^3/kg
    double pa = 50000;          // kg/m/s^2 (pascals)
    double aa = constant_state_soundspeed(rhoa,pa);

    //NOTE: ua most be supersonic, ub must be subsonic.

    printf("ua = %g\n",ua);
    printf("rhoa = %g\n",rhoa);
    printf("pa = %g\n",pa);
    printf("aa = %g\n\n",aa);

    //TODO: check STATE struct to see which variable has been provided?
    //double rhob = 1.225;
    //double taub = 1.0/rhob;
    double pb = 100000;

    //STATE sa = {ua, rhoa, pa};  //ahead state
    //STATE sb = {ub, rhob, pb};  //behind state

    //1.use Hugoniot function to compute pb OR rhob/taub (whichever not given)

    //double pb = behind_state_pressure(rhoa,pa,rhob);
    double taub = behind_state_specific_volume(rhoa,pa,pb);
    double rhob = 1.0/taub;
    double ab = constant_state_soundspeed(rhob,pb);

    //NOTE: M is mass flux in the shock-stationary frame, and is signed.

    //2.compute M from M^2 = (pr - pl)/(taul - taur)
    //2.compute M from M^2 = (pb - pa)/(taua - taub)

    double M = sign*std::sqrt((pb - pa)/(taua - taub));
    printf("M = %g\n\n",M);

    //3.set M = (pr - pl)/(ul - ur) and compute ur
    //3.set M = (pb - pa)/(ua - ub) and compute ub

    double ub = ua - (pb - pa)/M;

    printf("ub = %g\n",ub);
    printf("rhob = %g\n",rhob);
    printf("pb = %g\n",pb);
    printf("ab = %g\n\n",ab);

    //TODO: check prandtl relation, ua and ub have same sign for shock

    //4. check that entropy condition satisfied;
    //  i.e. rhoa < rhob && pa < pb

    bool entropy_condition = true;
    if (rhoa >= rhob)
    {
        entropy_condition = false;
        printf("ahead state density not less than behind state density\n");
    }
    if (pa >= pb)
    {
        entropy_condition = false;
        printf("ahead state pressure not less than behind state pressure\n");
    }
    if (!entropy_condition)
    {
        printf("entropy condition not satisfied\n");
    }

    double shock_speed = (rhoa*ua - rhob*ub)/(rhoa - rhob);
    double mach = (ua - shock_speed)/aa;

    printf("shock_speed = %g\n",shock_speed);
    printf("mach = %g\n",mach);

    //TODO: compute shock line and give solution for any (x,t) pair;
    //      locate region of point and compute solution.

    return 0;
}

