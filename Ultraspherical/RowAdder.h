//
//  RowAdder.h
//  Ultraspherical
//
//  Created by Sheehan Olver on 8/10/13.
//  Copyright (c) 2013 School of Mathematics and Statistics, The University of Sydney. All rights reserved.
//

#ifndef Ultraspherical_RowAdder_h
#define Ultraspherical_RowAdder_h


#include <iostream>
#include <vector>

using namespace std; 



class RowAdder
{
public:
    RowAdder();
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

#endif
