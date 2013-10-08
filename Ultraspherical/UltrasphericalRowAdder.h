

#ifndef Ultraspherical_Multiplication_Matrix_h
#define Ultraspherical_Multiplication_Matrix_h

#include <iostream>
#include "FilledBandedMatrix.h" 
#include "RowAdder.h"
#include <vector>

using namespace std; 


class ConversionRowAdder : public RowAdder
{
    unsigned int from,to;
    //TODO: Implement other conversion operators
    
public:
    ConversionRowAdder(unsigned int l, unsigned int m);
    
    virtual double getEntry(long,long);
    virtual long leftIndex(unsigned long);
    virtual long rightIndex(unsigned long);
};


class DerivativeRowAdder : public RowAdder
{
    unsigned int from,to;
public:
    
    DerivativeRowAdder(unsigned int l, unsigned int m);
    
    virtual double getEntry(long,long);
    virtual long leftIndex(unsigned long);
    virtual long rightIndex(unsigned long);
};



class ToeplitzRowAdder : public RowAdder
{
    //TODO: Not really Toeplitz
    double *a;
    unsigned long length;
    long index;
    
public:
    ToeplitzRowAdder(vector<double> *, long);
    ToeplitzRowAdder(double *, unsigned long, long);
    virtual double getEntry(long,long);
    virtual long leftIndex(unsigned long);
    virtual long rightIndex(unsigned long);
};

class HankelRowAdder : public RowAdder
{
    vector<double> *a;
    
public:
    HankelRowAdder(vector<double> *a);
    virtual double getEntry(long,long);
    virtual long leftIndex(unsigned long);
    virtual long rightIndex(unsigned long);
};


RowAdder *MultiplicationRowAdder(unsigned int, vector<double> *);


FilledBandedMatrix *DirichletD2ConvertMultiplicationMatrix(vector<double> *a);





#endif
