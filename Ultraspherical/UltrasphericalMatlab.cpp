//
//  Matlab.cpp
//  Ultraspherical
//
//  Created by Sheehan Olver on 22/07/2012.
//  Copyright (c) 2012 School of Mathematics and Statistics, The University of Sydney. All rights reserved.
//

//#include <iostream>

#include "mex.h"

// uncomment below and compile with mex main.cpp to generate MATLAB interface.

void mexFunction ( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    int argc; argc = 1; 
    const char * argv[1]; 
    //vector<double>* c; 
    
    //run main. 
    // main(argc,argv);
    
    //First entry is lower row, it's important!
//    FilledBandedOperator bnd(-2,&D2Dirichlet);
//    bnd.increaseSize();
//    bnd.increaseSize();
//    bnd.increaseSize();    
//    bnd.increaseSize();    
//    bnd.increaseSize(); 
//    bnd.increaseSize(); 
//    vector<double> b,c;
//    
//    b.push_back(1);
//    for(long i = 0; i < 50000; i++)
//        b.push_back(1);
//    
//    c = QRSolve(bnd,b);
    
    // return the coefficients of the solution. 
    const int len = (int) 2; 
    plhs[0] = mxCreateDoubleMatrix(len,1,mxREAL); // real coeffs vector. 
    double * y = mxGetPr(plhs[0]); 
    
    // assign coeff vector to output. 
    for( int i = 0; i<len; ++i)
        y[i] = nrhs;  
}
