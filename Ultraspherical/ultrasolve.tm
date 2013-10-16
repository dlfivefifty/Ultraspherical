
void ultraSolve P(( double *, long, double *, long, double *, long, double *, long));

:Begin:
:Function:       ultraSolve
:Pattern:        UltrasphericalSolve[i_List,j_List,bc_List,fn_List]
:Arguments:      { i, j, bc, fn }
:ArgumentTypes:  { RealList, RealList, RealList, RealList }
:ReturnType:     Manual
:End:

:Evaluate: UltrasphericalSolve::usage = "UltrasphericalSolve[{a0,...,an},{b0,...,bn},{bc0,bc1},{f0,...,fn}] solves the ODE u'' + (a0 T0(x) + ... + an T_n(x)) u' + (b0 T0(x) + ... + bn T_n(x)) u = f0 T0(x) + ... + fn T_n(x), u(-1) = bc0, u(1) = bc1."


void poissonSolve P(( double *, long, double *, long, double *, long, double *, long));

:Begin:
:Function:       poissonSolve
:Pattern:        PoissonSolve[i_List,j_List,bc_List,fn_List]
:Arguments:      { i, j, bc, fn }
:ArgumentTypes:  { RealList, RealList, RealList, RealList }
:ReturnType:     Manual
:End:

:Evaluate: PoissonSolve::usage = "UltrasphericalSolve[{a0,...,an},{b0,...,bn},{bc0,bc1},{f0,...,fn}] solves the ODE u'' + (a0 T0(x) + ... + an T_n(x)) u' + (b0 T0(x) + ... + bn T_n(x)) u = f0 T0(x) + ... + fn T_n(x), u(-1) = bc0, u(1) = bc1."
