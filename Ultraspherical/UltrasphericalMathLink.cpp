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
#include "MultiplicationMatrices.h"
//
////#include "mex.h"
#include "AdaptiveQR.h"

#include <vector>


extern "C" void ultraSolve( double *, long);



void ultraSolve( double *ain, long n, double *bin, long m)
{
    vector<double> a, b;
    
    for(int i = 0; i < n; i++)
        a.push_back(ain[i]);
    
    for (int i = 0; i < m; ++i) {
        b.push_back(bin[i]);
    }
    
    FilledBandedMatrix *drbnd = DirichletD2ConvertMultiplicationMatrix(&a,&b);
    drbnd->increaseSize();
    
    vector<double> f;
    
    f.push_back(1);
    f.push_back(0);
    f.push_back(0);
    
    vector<double> c = QRSolve(drbnd,f);   
    
    
    double cret[c.size()];
    
    for(int i = 0; i < c.size(); i++)
        cret[i] = c[i];
    
    
    
    MLPutReal64List(stdlink, cret, (int)c.size());
    
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
