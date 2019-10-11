#include "riemann_problem.h"


//LFS: moves to right with ul, ur > 0
//ul is ua and ur is ub

int main(int argc, char* argv[])
{
    //double gamma = 1.4;

    //Given behind and ahead states (sb and sa)
    //  use left and right for time being .. sl and sr
    
    double ul = 0.0;
    double rhol = 1.225;
    double pl = 100000;

    double ur = -1.0;
    double rhor = 1.1;
    double pr = 100;

    STATE sl = {ul, rhol, pl};
    STATE sr = {ur, rhor, pr};

    //given 4 inputs: ua, rhoa, pa and either rhob OR pb
    
    //1.use Hugoniot function to compute pb OR rhob (whichever not given)

    //note: M is not mach number, it is mass flux in the shock-stationary frame

    //2.compute M from M^2 = (p1 - p0)/(tau0 - tau1)
    //2.compute M from M^2 = (pr - pl)/(taul - taur)
    //2.compute M from M^2 = (pb - pa)/(taua - taub)

    //3.set M = (p1 - p0)/(u0 - u1) and compute u1
    //3.set M = (pr - pl)/(ul - ur) and compute ur
    //3.set M = (pb - pa)/(ua - ub) and compute ub

    //4. check that 
    return 0;
}
