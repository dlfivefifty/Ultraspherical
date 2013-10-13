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
    virtual double getEntry(unsigned long,unsigned  long);
    virtual long leftIndex(unsigned long);
    virtual long rightIndex(unsigned long);
    
    virtual Operator *operator+(Operator *);
    virtual Operator *operator*(double);
    
    void print();
};



class PlusOperator : public Operator
{
    vector<Operator *> *summands;
public:
    PlusOperator(Operator *rowAdd);
    ~PlusOperator();
    virtual double getEntry(unsigned  long,unsigned  long);
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
    virtual double getEntry(unsigned  long,unsigned  long);
    virtual long leftIndex(unsigned long);
    virtual long rightIndex(unsigned long);
};

class ConstantTimesOperator : public Operator
{
    double a;
    Operator *b;
public:
    ConstantTimesOperator(double a, Operator *bb);
    virtual double getEntry(unsigned long,unsigned long);
    virtual long leftIndex(unsigned long);
    virtual long rightIndex(unsigned long);
};


class ShiftOperator : public Operator
{
    Operator *adder;
    long shift;
public:
    ShiftOperator(Operator *, long);
    
    virtual double getEntry(unsigned long,unsigned long);
    virtual long leftIndex(unsigned long);
    virtual long rightIndex(unsigned long);
};

class SkipOperator : public Operator
{
    Operator *op;
    unsigned long row_start,row_skip,col_start,col_skip;
    
public:
    SkipOperator(Operator *,unsigned long,unsigned long,unsigned long,unsigned long);
    
    unsigned long shiftRow(unsigned long);
    unsigned long shiftColumn(unsigned long);
    
    virtual double getEntry(unsigned long,unsigned long);
    virtual long leftIndex(unsigned long);
    virtual long rightIndex(unsigned long);
};




#endif
