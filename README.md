Ultraspherical
==============

Implements ultraspherical spectral methods for solving linear ODEs. 
This method results in matrices which are almost banded and a solver
requires O(n) operations.


The package contains a C++, Matlab and Mathematica implementation.  All implementations are experimental.

The Matlab implementation requires chebfun:

        http://www2.maths.ox.ac.uk/chebfun/
        
The Mathematica implementation requires RHPackage:

        http://www.maths.usyd.edu.au/u/olver/projects/RHPackage.html
    
The C++ implementation can be accessed from Mathematica (without RHPackage).  See RunMathLink.nb.


Files:


ultraop                    		Matlab implementation folder
ultraop/examples        	Example usage of ultraop
INSTALL                 		Compilation and usage instructions for C++ implementation
RunMathLink.nb		Access C++ from Mathematica
init.nb/init.m			Mathematica implementation
Test                    		Unit tests


