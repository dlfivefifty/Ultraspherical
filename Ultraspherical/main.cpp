//
//  main.cpp
//  Ultraspherical
//
//  Created by Sheehan Olver on 04/06/2012.
//  Copyright (c) 2012 School of Mathematics and Statistics, The University of Sydney. All rights reserved.
//



#include "FilledBandedOperator.h" 
#include "UltrasphericalOperator.h"
#include "OlversAlgorithm.h"
#include "AdaptiveQR.h"
#include <time.h>
#include <math.h>


//void smallbandExample()
//{
//    cout << "SMALL BAND" <<endl;
//    
//    clock_t t1; clock_t t2;
//    
//    vector<double> b;
//    
//    for(long i = 0; i < 10; i++)
//        b.push_back(1);
//    
//    cout<<"{"<<endl;        
//    for(long n = 0; n < 41; n ++)
//    {
//        
//        vector<double> a;         
//        a.push_back(10); a.push_back(2); a.push_back(3);        
//        
//        FilledBandedOperator *drbnd = DirichletD2ConvertMultiplicationMatrix(NULL,&a);
//        
//        
//        //        drbnd.print();
//        
//        
//        drbnd->increaseSize();
//        
//        //        
//        for(long i = 0; i < n; i++)
//            b.push_back(1);
//        //    printvec(b);
//        
//        t1 = clock();                
//        
//        vector<double> c = QRSolve(drbnd,b);
//        
//        t2 = clock();
//        //            printvec(b);
//        float tottime = ((float)(t2-t1)/CLOCKS_PER_SEC);
//        cout<<"{"<<c.size()<<", "<<tottime<<"},"<<endl;  
//        
//        for(long i = 0; i < 50000; i++)
//            b.push_back(1);
//        
//        delete drbnd;
//    }
//    cout<<"}"<<endl;
//}


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


void smallbandExample()
{
    cout << "SMALL BAND" <<endl;
    
    clock_t t1; clock_t t2;
    
    long n = 40;
    double bin[n];
    
    for (int i = 0; i < n; ++i)
        bin[i] = 1.;
    
    toC2Series(bin, n);
    
    
    vector<double> b;
    
    b.push_back(1);
    b.push_back(1);
    
    for(long i = 0; i < n; i++)
        b.push_back(bin[i]);
    
    printvec(b);
    
    vector<double> a;
    a.push_back(10); a.push_back(2); a.push_back(3);
    
    FilledBandedOperator *drbnd = DirichletD2ConvertMultiplicationMatrix(NULL,&a);
    
    
    //        drbnd.print();
    
    
    drbnd->increaseSize();
    
    t1 = clock();
    
    vector<double> c = *QRSolve(drbnd,b);
    
    t2 = clock();
    //            printvec(b);
    float tottime = ((float)(t2-t1)/CLOCKS_PER_SEC);
    cout<<"{"<<c.size()<<", "<<tottime<<"},"<<endl;
    
    printvec(c);
    
    delete drbnd;
}

void cosbandExample()
{
    cout << "COSINE BAND" <<endl;    
    clock_t t1; clock_t t2;
    
    vector<double> b;
    
    for(long i = 0; i < 10; i++)
        b.push_back(1);
    
    cout<<"{"<<endl;        
    for(long n = 0; n < 41; n ++)
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
        
        FilledBandedOperator *drbnd = DirichletD2ConvertMultiplicationMatrix(NULL, &a);
        drbnd->increaseSize();
        
        //        

        //    printvec(b);
        
        t1 = clock();                
        
        vector<double> c = *QRSolve(drbnd,b);
        
        t2 = clock();
        //            printvec(b);
        float tottime = ((float)(t2-t1)/CLOCKS_PER_SEC);
        cout<<"{"<<c.size()<<", "<<tottime<<"},"<<endl;    
        
        
        for(long i = 0; i < 50000; i++)
            b.push_back(1);
        
        delete drbnd;
    }
    cout<<"}"<<endl;
}


void airyExample()
{
    int numex = 14;
    double left[] = {0.5355608832923521, 0.1271728034584682, 0.351913571284756, 
        0.04024123848644319, -0.2607345878897477, -0.19426241417077472, 
        0.1767533932395529, -0.12078802581383595, 0.10178242352993307, 
        0.05597189577301992, 0.02333382924837296,-0.027905156151965354, 0.02705738360464258, 0.044775817580242405, 
            0.01829406779298808, -0.013152978737498166, -0.02249517692130777, 
            0.020321506196299167, -0.0021912611413430574, 
            -0.013261784699976999, 0.005242247212492905}  ;
    double right[] = {0.13529241631288141, 0.027515172874532375,
        0.00024221691467447962, 1.1047532552898686E-10, 
        1.4576297592861973E-30, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    for (int k = 0; k < numex; k++) {        
        clock_t t1; clock_t t2;
        vector<double> a;
        a.push_back(0);
        a.push_back(-pow(10, k));
        
        vector<double> b;
        b.push_back(left[k]);
        b.push_back(right[k]);
        
        FilledBandedOperator *drbnd = DirichletD2ConvertMultiplicationMatrix(NULL,&a);
        drbnd->increaseSize();
        
        t1 = clock();                
        
        vector<double> c = *QRSolve(drbnd,b);
        
        t2 = clock();
        //    printvec(c);
        float tottime = ((float)(t2-t1)/CLOCKS_PER_SEC);
        cout<<"{"<<c.size()<<", "<<tottime<<"},"<<endl;
        
        delete drbnd;
    }
}



void olversAlgorithm()
{
    vector<double> f;
    
    f.push_back(0.0625);
    f.push_back(0.0625);

    vector<vector<double> *> *uret = poisson(&f);
    
    int i = 0;
    for (vector<double> * ui  : *uret) {
        cout <<endl << "u[" << i <<"]: "<< endl;
        printvec(*ui);
        ++i;
    }
//
    
    
//    for (unsigned long i = 0; i < OpLength; ++i) {
////        printvec(*r[i]);
//        delete A[i];
//    }
}



int main(int argc, const char * argv[])
{
    
    if (argc > 1) {
        switch (argv[1][0]) {
            case 'c':
                cosbandExample();
                return 0;
            case 's':
                smallbandExample();
                return 0;
            case 'a':
                airyExample();
                return 0;
            case 'o':
                olversAlgorithm();
                return 0;
            default:
                break;
        }
    }
    
    if(argc < 4) {
        cout <<endl<< "Usage: "<<endl<<endl<<"\t\tultraspherical cm cp a0 a1 ... an"<<endl<<endl<<"solves the ODE "<<endl<<endl<<"\t\tu'' + (a0 T_0(x) + ... + an T_n(x)) u = 0, u(-1) = cm, u(1) = cp"<<endl<<endl;
        
        return 0;
    }
    
    
    vector<double> args;
    vector<double> b;
    double a;
    sscanf(argv[1],"%lf",&a);    
    b.push_back(a);
    sscanf(argv[2],"%lf",&a);    
    b.push_back(a);          
    b.push_back(0);
    
    for (int i = 3; i < argc; i++) {
        sscanf(argv[i],"%lf",&a);
        args.push_back(a);
    }
    

    
    FilledBandedOperator *drbnd = DirichletD2ConvertMultiplicationMatrix(&args, &args);

    
    drbnd->increaseSize();
//    drbnd->print();
    
    
    
    
    vector<double> c = *QRSolve(drbnd,b);
    
    printvec(c);
    
//    string args(argv[0]);
    
    

    

    
    delete drbnd;
    
    
//    cosbandExample();


    
    
    return 0;
}



