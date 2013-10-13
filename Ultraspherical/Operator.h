//
//  Operator.h
//  Ultraspherical
//
//  Created by Sheehan Olver on 8/10/13.
//  Copyright (c) 2013 School of Mathematics and Statistics, The University of Sydney. All rights reserved.
//

#ifndef Ultraspherical_Operator_h
#define Ultraspherical_Operator_h


#include <iostream>
#include <vector>

using namespace std; 



class Operator
{
public:
    Operator();
    virtual double getEntry(long,long);
    virtual long leftIndex(unsigned long);
    virtual long rightIndex(unsigned long);
    
    void print();
};



class PlusOperator : public Operator
{
    vector<Operator *> *summands;
public:
    PlusOperator(Operator *rowAdd);
    virtual double getEntry(long,long);
    virtual long leftIndex(unsigned long);
    virtual long rightIndex(unsigned long);
    void push_back(Operator *rowAdd);
    
};


class TimesOperator : public Operator
{
    Operator *a;
    Operator *b;
public:
    TimesOperator(Operator *aa, Operator *bb);
    virtual double getEntry(long,long);
    virtual long leftIndex(unsigned long);
    virtual long rightIndex(unsigned long);
};

class DoubleTimesOperator : public Operator
{
    double a;
    Operator *b;
public:
    DoubleTimesOperator(double a, Operator *bb);
    virtual double getEntry(long,long);
    virtual long leftIndex(unsigned long);
    virtual long rightIndex(unsigned long);
};


class ShiftOperator : public Operator
{
    Operator *adder;
    long shift;
public:
    ShiftOperator(Operator *, long);
    
    virtual double getEntry(long,long);
    virtual long leftIndex(unsigned long);
    virtual long rightIndex(unsigned long);
};




#endif
