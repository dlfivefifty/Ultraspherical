/* To launch this program from within Mathematica use:
 *   In[1]:= link = Install["addtwo"]
 *
 * Or, launch this program from a shell and establish a
 * peer-to-peer connection.  When given the prompt Create Link:
 * type a port name. ( On Unix platforms, a port name is a
 * number less than 65536.  On Mac or Windows platforms,
 * it's an arbitrary word.)
 * Then, from within Mathematica use:
 *   In[1]:= link = Install["portname", LinkMode->Connect]
 */

#include "mathlink.h"

#include "FilledBandedMatrix.h" // stop multiple definition of RowFiller class. 
#include "UltrasphericalRowAdder.h"
//
////#include "mex.h"
#include "AdaptiveQR.h"

#include <vector>


extern "C" void ultraSolve( double *, long, double *, long, double *, long, double *, long);


void toUSeries(double *fin, long fn)
{
    switch (fn) {
        case 0:
        case 1:
            break;
            
        case 2:
            fin[1] = .5*fin[1];
            break;
            
        default:
            fin[0] = fin[0] - .5*fin[2];
            
            for (int i = 1; i < fn-2; ++i) {
                fin[i] = .5*fin[i] - .5*fin[i+2];
            }
            
            fin[fn-2] = .5*fin[fn-2];
            fin[fn-1] = .5*fin[fn-1];
            break;
    }

}

void toC2Series(double *fin, long fn)
{
    toUSeries(fin, fn);
    
    switch (fn) {
        case 0:
        case 1:
            break;
            
        case 2:
            fin[1] = .5*fin[1];
            break;
            
        default:
            for (int i = 0; i < fn-2; ++i) {
                fin[i] = 1/((double) i+1)*fin[i] -  1/((double) i+3)*fin[i+2];
            }
            
            fin[fn-2] = 1/((double) fn-1)*fin[fn-2];
            fin[fn-1] = 1/((double) fn)*fin[fn-1];
            break;
    }
    
}


void ultraSolve( double *ain, long n, double *bin, long m, double *bc, long k, double *fin, long fn)
{
    vector<double> a, b;
    
    if (k != 2)
        return;
    
    for(int i = 0; i < n; i++)
        a.push_back(ain[i]);
    
    for (int i = 0; i < m; ++i) {
        b.push_back(bin[i]);
    }
    
    FilledBandedMatrix *drbnd = DirichletD2ConvertMultiplicationMatrix(&a,&b);
    drbnd->increaseSize();
    
    vector<double> f;
    
    f.push_back(bc[0]);
    f.push_back(bc[1]);

    toC2Series(fin,fn);
    
    
    
    for (long i = 0; i < fn; ++i) {
        f.push_back(fin[i]);
    }

    
    vector<double> c = QRSolve(drbnd,f);   
    
    
    double cret[c.size()];
    
    for(long i = 0; i < c.size(); i++)
        cret[i] = c[i];
    
    
    
    MLPutReal64List(stdlink, cret, (int)c.size());
    
    delete drbnd;
    
    return ;
//    return c[0];
}



#if WINDOWS_MATHLINK

#if __BORLANDC__
#pragma argsused
#endif

int PASCAL WinMain( HINSTANCE hinstCurrent, HINSTANCE hinstPrevious, LPSTR lpszCmdLine, int nCmdShow)
{
	char  buff[512];
	char FAR * buff_start = buff;
	char FAR * argv[32];
	char FAR * FAR * argv_end = argv + 32;
    
	hinstPrevious = hinstPrevious; /* suppress warning */
    
	if( !MLInitializeIcon( hinstCurrent, nCmdShow)) return 1;
	MLScanString( argv, &argv_end, &lpszCmdLine, &buff_start);
	return MLMain( (int)(argv_end - argv), argv);
}

#else

int main(int argc, char* argv[])
{
	return MLMain(argc, argv);
}

#endif
