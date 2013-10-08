

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
    FilledRow *createRow(unsigned long);
    virtual double getEntry(long,long);
    virtual long leftIndex(unsigned long);
    virtual long rightIndex(unsigned long);
};



class ConversionRowAdder : public RowAdder
{
//TODO: Implement other conversion operators
    
public:
    virtual double getEntry(long,long);
    virtual long leftIndex(unsigned long);
    virtual long rightIndex(unsigned long);
};


class PlusRowAdder : public RowAdder
{
    vector<RowAdder *> *summands;
public:
    PlusRowAdder(RowAdder *rowAdd);
    virtual double getEntry(long,long);
    virtual long leftIndex(unsigned long);
    virtual long rightIndex(unsigned long);
    void push_back(RowAdder *rowAdd);
    
};


class TimesRowAdder : public RowAdder
{
    vector<RowAdder *> *summands;
public:
    TimesRowAdder(RowAdder *rowAdd);
    virtual double getEntry(long,long);
    virtual long leftIndex(unsigned long);
    virtual long rightIndex(unsigned long);
    void push_back(RowAdder *rowAdd);
};


class DerivativeRowAdder : public RowAdder
{
    //TODO: Implement other derivative operators
public:
    virtual double getEntry(long,long);
    virtual long leftIndex(unsigned long);
    virtual long rightIndex(unsigned long);
};

class ShiftRowAdder : public RowAdder
{
    RowAdder *adder;
    long shift;
public:
    ShiftRowAdder(RowAdder *, long);
    
    virtual double getEntry(long,long);
    virtual long leftIndex(unsigned long);
    virtual long rightIndex(unsigned long);
};

class DirichletD2ConvertMultiplicationMatrix : public FilledBandedMatrix
{
    ShiftRowAdder *adder;
    
public:
    DirichletD2ConvertMultiplicationMatrix(DirichletD2ConvertMultiplicationMatrix& other);
    DirichletD2ConvertMultiplicationMatrix(vector<double> a);
    virtual FilledRow *createRow(unsigned long);
};






#endif
