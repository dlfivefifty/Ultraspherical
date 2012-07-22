
double ultraSolve P(( double *, long));

:Begin:
:Function:       ultraSolve
:Pattern:        UltrasphericalSolve[i_List]
:Arguments:      { i }
:ArgumentTypes:  { RealList }
:ReturnType:     Real
:End:

:Evaluate: UltrasphericalSolve::usage = "UltrasphericalSolve[{a0,...,an}] solves the ODE u'' + (a0 T0(x) + ... + an Tn(x)) u = 0."
