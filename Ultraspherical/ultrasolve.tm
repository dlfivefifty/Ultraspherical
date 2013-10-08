
void ultraSolve P(( double *, long, double *, long));

:Begin:
:Function:       ultraSolve
:Pattern:        UltrasphericalSolve[i_List,j_List]
:Arguments:      { i, j }
:ArgumentTypes:  { RealList, RealList }
:ReturnType:     Manual
:End:

:Evaluate: UltrasphericalSolve::usage = "UltrasphericalSolve[{a0,...,an}] solves the ODE u'' + (a0 T0(x) + ... + an Tn(x)) u = 0, u(-1) = 1, u(1) = 0."
