//
//  main.cpp
//  Ultraspherical
//
//  Created by Sheehan Olver on 04/06/2012.
//  Copyright (c) 2012 School of Mathematics and Statistics, The University of Sydney. All rights reserved.
//


//#include "FilledBandedMatrix.h"
#include "FilledBandedMatrix.h" // stop multiple definition of RowFiller class. 
#include "MultiplicationMatrices.h"

//#include "mex.h"
#include "AdaptiveQR.h"
#include <time.h>



void smallbandExample()
{
    clock_t t1; clock_t t2;
    
    vector<double> b;
    
    b.push_back(1);    
    
    cout<<"{"<<endl;        
    for(long n = 10; n < 11; n += 100)
    {
        
        vector<double> a;         
        a.push_back(10); a.push_back(2); a.push_back(3);        
        
        DirichletD2ConvertMultiplicationMatrix drbnd(a);
        
        
        //        drbnd.print();
        
        
        drbnd.increaseSize();
        
        //        
        for(long i = 0; i < n; i++)
            b.push_back(1);
        //    printvec(b);
        
        t1 = clock();                
        
        vector<double> c = QRSolve(&drbnd,b);            
        
        t2 = clock();
        //            printvec(b);
        float tottime = ((float)(t2-t1)/CLOCKS_PER_SEC);
        cout<<"{"<<c.size()<<", "<<tottime<<"},"<<endl;    
    }
    cout<<"}"<<endl;              
}

void cosbandExample()
{
    clock_t t1; clock_t t2;
    
    vector<double> b;
    
    for(long i = 0; i < 10; i++)
        b.push_back(1);
    
    cout<<"{"<<endl;        
    for(long n = 0; n < 10; n ++)
    {
        
        vector<double> a;         
        a.push_back(0.7651976865579668);
        a.push_back(0.);
        a.push_back(-0.22980696986380098);
        a.push_back(0.);
        a.push_back(0.004953277928219901);
        a.push_back(0.);
        a.push_back(-0.00004187667600477857);
        a.push_back(0.);
        a.push_back(1.8844688343861268e-7);
        a.push_back(0.);
        a.push_back(-5.261244076090109e-10);
        a.push_back(0.);a.push_back(9.999074386172517e-13);  
        
        DirichletD2ConvertMultiplicationMatrix drbnd(a);             
        drbnd.increaseSize();
        
        //        

        //    printvec(b);
        
        t1 = clock();                
        
        vector<double> c = QRSolve(&drbnd,b);            
        
        t2 = clock();
        //            printvec(b);
        float tottime = ((float)(t2-t1)/CLOCKS_PER_SEC);
        cout<<"{"<<c.size()<<", "<<tottime<<"},"<<endl;    
        
        
        for(long i = 0; i < 50000; i++)
            b.push_back(1);        
    }
    cout<<"}"<<endl;              
}


void airyExample()
{
//    double e = 1/10.;
//    double left = 0.1271728034584682;
//    double right = 0.027515172874532375;
    
//    double e = 1/100.;
//    double left = 0.351913571284756;
//    double right = 0.00024221691467447962;    

//    int numex = 10;
    double e = 1./1000000000;
    double left = 0.05597189577301992;
    double right = 0;
    
    
    clock_t t1; clock_t t2;
    vector<double> a;
    a.push_back(0);
    a.push_back(-1/e);
    
    vector<double> b;
    b.push_back(left);
    b.push_back(right);
    
    DirichletD2ConvertMultiplicationMatrix drbnd(a);             
    drbnd.increaseSize();

    t1 = clock();                
    
    vector<double> c = QRSolve(&drbnd,b);            
    
    t2 = clock();
//    printvec(c);
    float tottime = ((float)(t2-t1)/CLOCKS_PER_SEC);
    cout<<"{"<<c.size()<<", "<<tottime<<"},"<<endl;       
}



int main(int argc, const char * argv[])
{    

    airyExample();


    
    return 0;
}

// uncomment below and compile with mex main.cpp to generate MATLAB interface.
/*
void mexFunction ( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

int argc; argc = 1; 
const char * argv[1]; 
//vector<double>* c; 

//run main. 
main(argc,argv);

    //First entry is lower row, it's important!
    FilledBandedMatrix bnd(-2,&D2Dirichlet);
        bnd.increaseSize();
        bnd.increaseSize();
        bnd.increaseSize();    
        bnd.increaseSize();    
        bnd.increaseSize(); 
        bnd.increaseSize(); 
    vector<double> b,c;

    b.push_back(1);
    for(long i = 0; i < 50000; i++)
        b.push_back(1);
    
    c = QRSolve(bnd,b);

// return the coefficients of the solution. 
const int len = (int)c.size(); 
plhs[0] = mxCreateDoubleMatrix(len,1,mxREAL); // real coeffs vector. 
double * y = mxGetPr(plhs[0]); 

// assign coeff vector to output. 
for( int i = 0; i<len; ++i)
    y[i] = c[i];  
}
*/

