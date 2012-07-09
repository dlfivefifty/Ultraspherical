//
//  main.cpp
//  Ultraspherical
//
//  Created by Sheehan Olver on 04/06/2012.
//  Copyright (c) 2012 School of Mathematics and Statistics, The University of Sydney. All rights reserved.
//


//#include "FilledBandedMatrix.h"
#include "FilledBandedMatrix_extras.h" // stop multiple definition of RowFiller class. 
#include "MultiplicationMatrices.h"

//#include "mex.h"
#include "AdaptiveQR.h"
#include <time.h>


//FilledRow defaultRowCreator(int k)
//{
//    
//    double w = 1;
//    
//    vector<RowFiller> fillers;
//    fillers.push_back(RowFiller::leftDirichlet());
//    fillers.push_back(RowFiller::rightDirichlet());    
//    
//    switch(k)
//    {
//        case 0:
//        {
//            FilledRow newrow(0,RowFiller::dirichlet(1,0));
//            newrow.increaseSize();//must have an element in first row to use Givens                        
//            return newrow;
//        }
//        case 1:
//        {
//            FilledRow newrow(0,RowFiller::dirichlet(0,1));
//            newrow.increaseSize();//must have an element in first row to delete            
//            return newrow;
//        }            
//        case 2:
//        {
//            FilledRow newrow(0,RowFiller::dirichlet(0,0),3);               
//            newrow.setEntry(0,w);
//            newrow.setEntry(1,1);
//            newrow.setEntry(2,-.5*w);
//        
//            return newrow;
//        }
//        default:
//        {
//            FilledRow newrow(k-2,RowFiller::dirichlet(0,0),3);       
//            
//            newrow.setEntry(k-2,.5*w);
//            newrow.setEntry(k-1,k-1);
//            newrow.setEntry(k,-.5*w);            
//            
//            return newrow;
//        }
//    }    
//}


FilledRow D2Dirichlet(int k)
{
    
//    double w = 1;
    
    switch(k)
    {
        case 0://LeftDir
        {
            FilledRow newrow(0,RowFiller::dirichlet(1,0));
            newrow.increaseSize();//must have an element in first row to use Givens                        
            return newrow;
        }
        case 1://RightDir
        {
            FilledRow newrow(0,RowFiller::dirichlet(0,1));
            newrow.increaseSize();//must have an element in first row to delete            
            return newrow;
        }     
        case 2:
        {
            FilledRow newrow(0,RowFiller::dirichlet(0,0));
            newrow.push_back(1);
            newrow.push_back(0);            
            newrow.push_back(4.-2./3.);            
            newrow.push_back(0);            
            newrow.push_back(1./6);                        
            return newrow;    
        }
        default:
        {
            FilledRow newrow(k-2,RowFiller::dirichlet(0,0));
            newrow.push_back(1./(2.*(k-1.)));
            newrow.push_back(0);            
            newrow.push_back(2.*k-k/(2.*(k-1)+(k-1)*(k-1)));            
            newrow.push_back(0);            
            newrow.push_back(1./(2.*(k+1)));                   
            return newrow;  
        }
    }    
}


int main(int argc, const char * argv[])
{    
    //First entry is lower row, it's important!
    FilledBandedMatrix bnd(-2,&D2Dirichlet);
    clock_t t1; clock_t t2;

    
    
    
    
        bnd.increaseSize();
        bnd.increaseSize();
        bnd.increaseSize();    
        bnd.increaseSize();    
        bnd.increaseSize(); 
        bnd.increaseSize();     
    
    
    cout<<"B: "<<endl;
    bnd.print();
    //    
//        cout<<"GB: "<<endl;    
//        bnd.applyGivens(0,1,0,1);
//        bnd.print();
//    
    
    //    
    //    
    vector<double> b;

    b.push_back(1);
    

    vector<double> a; 
//    double omega = -10;
//    a.push_back(omega); a.push_back(2*omega); a.push_back(3*omega); //a.push_back(3);  a.push_back(4);  a.push_back(5); a.push_back(6);a.push_back(7);
  

    //  Cosine
//    a.push_back(0.7651976865579668);
//    a.push_back(0.);
//    a.push_back(-0.22980696986380098);
//    a.push_back(0.);
//    a.push_back(0.004953277928219901);
//    a.push_back(0.);
//    a.push_back(-0.00004187667600477857);
//    a.push_back(0.);
//    a.push_back(1.8844688343861268e-7);
//    a.push_back(0.);
//    a.push_back(-5.261244076090109e-10);
//    a.push_back(0.);a.push_back(9.999074386172517e-13);
//    cout <<"toep:"<<endl;
    
    //    hankelOperator(a).print();
//    ConvertToeplitzMatrix(a).print();
//    cout <<"hankel:"<<endl;    
//    ConvertHankelMatrix(a).print();    
//    cout <<"mult:"<<endl;        
    
    
//    drbnd.print();        
    
//    cout<<"lower:"<<drbnd.lower()<<endl;
    DirichletD2ConvertMultiplicationMatrix(a).print();


//    
        cout<<"{"<<endl;        
    for(long n = 0; n < 10; n ++)
    {

        
        DirichletD2ConvertMultiplicationMatrix drbnd(a);
        
        
        drbnd.increaseSize();

        
        for(long i = 0; i < 100*n*n; i++)
            b.push_back(1);
        //    printvec(b);
        
        t1 = clock();                
        
        vector<double> c = QRSolve(&drbnd,b);            
        
        t2 = clock();
        //    printvec(b);
        //    printvec(c);    
        float tottime = ((float)(t2-t1)/CLOCKS_PER_SEC);
        cout<<"{"<<c.size()<<", "<<tottime<<"},"<<endl;    
    }
        cout<<"}"<<endl;            
        
    cout << b.size() <<endl;
    
//    cout<<"time/size: "<<tottime*300000/b.size()<<endl;    


    
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

