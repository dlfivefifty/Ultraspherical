

#ifndef Ultraspherical_Multiplication_Matrix_h
#define Ultraspherical_Multiplication_Matrix_h

#include <iostream>
#include "FilledBandedMatrix.h" 
#include "RowAdder.h"
#include <vector>

using namespace std; 



class ToeplitzRowAdder : public RowAdder
{
    //TODO: Not really Toeplitz
    vector<double> *a;
    
public:
    ToeplitzRowAdder(vector<double> *a);
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


RowAdder *MultiplicationRowAdder(vector<double> *);


FilledBandedMatrix *DirichletD2ConvertMultiplicationMatrix(vector<double> *a);





#endif
