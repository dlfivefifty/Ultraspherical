

#ifndef Ultraspherical_Multiplication_Matrix_h
#define Ultraspherical_Multiplication_Matrix_h

#include <iostream>
#include "FilledBandedMatrix.h" 
#include <vector>

using namespace std; 

//class ToeplitzMatrix : public FilledBandedMatrix
//{
//    vector<double> *rowEntries;
//public:
//    ToeplitzMatrix(vector<double> a);
//    virtual FilledRow createRow(int);
//};
//
//class MultiplicationMatrix : public FilledBandedMatrix
//{
//    vector<double> *rowEntries;
//public:
//    MultiplicationMatrix(vector<double> a);
//    virtual FilledRow createRow(int);
//};
//
//class ConvertToeplitzMatrix : public FilledBandedMatrix
//{
//    vector<double> *rowEntries;
//public:
//    ConvertToeplitzMatrix(vector<double> a);
//    virtual FilledRow createRow(int);
//};
//
//
//class ConvertHankelMatrix : public FilledBandedMatrix
//{
//public:
//    ConvertHankelMatrix(vector<double> a);
//    virtual FilledRow createRow(int);
//};
//
//class ConvertMultiplicationMatrix : public FilledBandedMatrix
//{
//    vector<double> *rowEntries;
//public:
//    ConvertMultiplicationMatrix(vector<double> a);
//    virtual FilledRow createRow(int);
//};

class RowAdder
{
public:
    RowAdder();
    virtual FilledRow *createRow(unsigned long);
};



class DirichletD2ConvertMultiplicationRowAdder : public RowAdder
{
    vector<double> *rowEntries;
    
public:
    DirichletD2ConvertMultiplicationRowAdder(vector<double> *a);
    virtual FilledRow *createRow(unsigned long);
};


class PlusRowAdder : public RowAdder
{
    vector<RowAdder *> *summands;
public:
    PlusRowAdder(RowAdder *rowAdd);
    virtual FilledRow *createRow(unsigned long);
    void push_back(RowAdder *rowAdd);
    
};

//only second currently
class DerivativeRowAdder : public RowAdder
{
    virtual FilledRow *createRow(unsigned long);
};

class DirichletD2ConvertMultiplicationMatrix : public FilledBandedMatrix
{
    PlusRowAdder *adder;
    
public:
    DirichletD2ConvertMultiplicationMatrix(DirichletD2ConvertMultiplicationMatrix& other);
    DirichletD2ConvertMultiplicationMatrix(vector<double> a);
    virtual FilledRow *createRow(unsigned long);
};






#endif
