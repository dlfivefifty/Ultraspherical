

#ifndef Ultraspherical_Multiplication_Matrix_h
#define Ultraspherical_Multiplication_Matrix_h

#include <iostream>
#include "FilledBandedOperator.h" 
#include "Operator.h"
#include <vector>

using namespace std; 


class ConversionOperator : public Operator
{
    unsigned int from,to;
    //TODO: Implement other conversion operators
    
public:
    ConversionOperator(unsigned int l, unsigned int m);
    
    virtual double getEntry(unsigned long,unsigned long);
    virtual long leftIndex(unsigned long);
    virtual long rightIndex(unsigned long);
};


class DerivativeOperator : public Operator
{
    unsigned int from,to;
public:
    
    DerivativeOperator(unsigned int l, unsigned int m);
    
    virtual double getEntry(unsigned long,unsigned long);
    virtual long leftIndex(unsigned long);
    virtual long rightIndex(unsigned long);
};



class ToeplitzOperator : public Operator
{
    //TODO: Not really Toeplitz
    double *a;
    unsigned long length;
    long index;
    
public:
    ToeplitzOperator(vector<double> *, long);
    ToeplitzOperator(double *, unsigned long, long);
    virtual double getEntry(unsigned long,unsigned long);
    virtual long leftIndex(unsigned long);
    virtual long rightIndex(unsigned long);
};

class HankelOperator : public Operator
{
    double *a;
    unsigned long length;
    
public:
    HankelOperator(vector<double> *a);
    HankelOperator(double *, unsigned long);
    virtual double getEntry(unsigned long,unsigned long);
    virtual long leftIndex(unsigned long);
    virtual long rightIndex(unsigned long);
};


Operator *MultiplicationOperator(unsigned int, vector<double> *);


FilledBandedOperator *DirichletD2ConvertMultiplicationMatrix(vector<double> *,vector<double> *);





#endif
