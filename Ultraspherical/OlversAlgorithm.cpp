//
//  OlversAlgorithm.cpp
//  Ultraspherical
//
//  Created by Sheehan Olver on 10/10/13.
//  Copyright (c) 2013 School of Mathematics and Statistics, The University of Sydney. All rights reserved.
//

#include "OlversAlgorithm.h"

#include "AdaptiveQR.h"

#include <math.h>



Operator *ConstantOperator(double c)
{
    double *one = (double *)malloc(sizeof(double));
    
    *one = c;
    
    return new ToeplitzOperator(one, 1,0);
}




RSOperator::RSOperator()
{
}

double RSOperator::getEntry(unsigned long row, unsigned long col)
{
    double c = (double) col + 2;
    double r = 1/((double)4 + 2*row);
    
    if (row == col + 2)
        return r/(2.*(c+1));
    else if (row +2 == col + 2)
        return r*(-1/(2*(c+1))-1/(2*(c-1)));
    else if (row + 4 == col + 2)
        return r/(2*(c-1));
    
                

    return 0;
}


long RSOperator::leftIndex(unsigned long row)
{
    if (row < 2) {
        return 0;
    } else
        return row - 2;
}

long RSOperator::rightIndex(unsigned long row)
{
    return row + 2;
}


Operator *RSListOperator(double *a, double *s, unsigned long n)
{
    if(n < 1)
        return NULL;
        
    
    PlusOperator *pl = new PlusOperator(new ConstantTimesOperator(a[0], new RSOperator()));
    pl->push_back(ConstantOperator(s[0]));
    
    
    
    return pl;
}


vector<double> *vectorTimes(vector<double> *a, double c)
{
    vector<double> *ret = new vector<double>;
    for (double en : *a)
    {
        ret->push_back(c*en);
    }
    
    return ret;
}



vector<double> *vectorPlus(vector<double> *a, vector<double> *b)
{
    vector<double> *ret = new vector<double>;
    
    unsigned long n = a->size(),m = b->size();
    
    for (unsigned long i = 0; i < max(n,m); ++i) {
        double en = 0;
        if(i < n)
            en += (*a)[i];
        
        if (i < m)
            en += (*b)[i];
        
        ret->push_back(en);
    }
    
    return ret;
}


double norm(vector<double> *a)
{
    double ret = 0;
    
    for(double i : *a)
        ret += i*i;
    
    return sqrt(ret);
}


vector<vector<double> *> *poisson(vector<double> *f)
{
    
    
    Operator *S = new SavedOperator(new SkipOperator(new RSOperator(), 0, 2, 0, 2));
    Operator *I = new SavedOperator(ConstantOperator(1));
    
#define beta(k) S->getEntry(k,k+1)
#define gamma(k) S->getEntry(k,k-1)
    
    
    vector<vector<double> *> r;
    
    
    r.push_back(vectorTimes(f, gamma(1)));
    r.push_back(vectorTimes(r[0], -1));
    
    
    vector<Operator *> A, B;
    
    A.push_back((*S) + (*I)*S->getEntry(0,0));
    A.push_back((*S) + (*I)*S->getEntry(1,1));
    B.push_back(A[0]);
    B.push_back(*((*B[0])*A[1]) + (*I)*(-beta(0)*gamma(1)));
    
#define dabs(d) d < 0? -d : d
    
    for (unsigned long i = 2;  norm(r[i-2])/dabs(B[i-2]->getEntry(0,0)) > 1.0E-16; ++i) {
        r[i - 1] =  vectorTimes(r[i-1], gamma(i));
        r.push_back(vectorTimes(r[i-1], -1));
        A.push_back((*S) + (*I)*S->getEntry(i,i));
        B.push_back(*((*B[i-1])*A[i]) + (*B[i-2])*(-beta(i-1)*gamma(i)));
    }
    
    const long OpLength = r.size();
    
    

    
    
    
   
    
    
    
#define Rsup(i) ((*B[i-1])*(beta(i)*gamma(i+1)))
#define Rdiag(i) (((*B[i])*gamma(i+1)))

    
    
    //Block back substitution
    vector<double> *u[OpLength - 1];
    //
    //
    int i = (int)OpLength - 2;
    u[i] = QRSolve(new FilledBandedOperator(-1-i, Rdiag(i)), *r[i]);
    //
    for (i = (int)OpLength - 3; i >= 1; --i) {
        r[i] = vectorPlus(r[i],vectorTimes((*Rsup(i))*(u[i+1]),-1));
//        cout<< "r:"<<endl;
//        printvec(*r[i]);
        u[i] = QRSolve(new FilledBandedOperator(-1-i, Rdiag(i)), *r[i]);
//        cout<<endl<< "u:"<<endl;
//        printvec(*u[i]);
    }
    
    i = 0;
    r[i] = vectorPlus(r[i],vectorTimes((u[i+1]),-beta(i)*gamma(i+1)));
    u[i] = QRSolve(new FilledBandedOperator(-1-i, Rdiag(i)), *r[i]);
    //

    
    vector<vector<double> *> *uret = new vector<vector<double> *>;
    
    for(long i = 0; i < OpLength - 1; ++i)
        uret->push_back(u[i]);
    
    return uret;
    
}


