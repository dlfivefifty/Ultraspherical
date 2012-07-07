

#ifndef Ultraspherical_Multiplication_Matrix_h
#define Ultraspherical_Multiplication_Matrix_h

#include <iostream>
#include "FilledBandedMatrix_extras.h" 
#include <vector>

using namespace std; 

class ToeplitzMatrix : public FilledBandedMatrix
{
    vector<double> rowEntries;
public:
    ToeplitzMatrix(vector<double> a);
    virtual FilledRow createRow(int);
};

class MultiplicationMatrix : public FilledBandedMatrix
{
    vector<double> rowEntries;
public:
    MultiplicationMatrix(vector<double> a);
    virtual FilledRow createRow(int);
};

class ConvertToeplitzMatrix : public FilledBandedMatrix
{
    vector<double> rowEntries;
public:
    ConvertToeplitzMatrix(vector<double> a);
    virtual FilledRow createRow(int);
};


class ConvertHankelMatrix : public FilledBandedMatrix
{
public:
    ConvertHankelMatrix(vector<double> a);
    virtual FilledRow createRow(int);
};

class ConvertMultiplicationMatrix : public FilledBandedMatrix
{
    vector<double> rowEntries;
public:
    ConvertMultiplicationMatrix(vector<double> a);
    virtual FilledRow createRow(int);
};

class DirichletD2ConvertMultiplicationMatrix : public FilledBandedMatrix
{
    vector<double> rowEntries;
public:
    DirichletD2ConvertMultiplicationMatrix(DirichletD2ConvertMultiplicationMatrix& other);
    DirichletD2ConvertMultiplicationMatrix(vector<double> a);
    virtual FilledRow createRow(int);
};




vector<double> ConvertMult(vector<double>,int,int);// S_0 * M_0[a]

vector<double> ConvertHankel(vector<double>,int); //used for convertmult

FilledBandedMatrix hankelOperator(vector<double>);

vector<double> HankelDiff(vector<double>, int);
//vector<double> MultDiff(vector<double>, int);




#endif
