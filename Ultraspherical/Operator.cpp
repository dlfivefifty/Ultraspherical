//
//  Operator.c
//  Ultraspherical
//
//  Created by Sheehan Olver on 8/10/13.
//  Copyright (c) 2013 School of Mathematics and Statistics, The University of Sydney. All rights reserved.
//

#include <stdio.h>

#include "Operator.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>



Operator::Operator()
{
    
}

double Operator::getEntry(unsigned long k, unsigned long j)
{
    cout << "getEntry NOT DEFINED!!"<<endl;
    return -1;
}

long Operator::leftIndex(unsigned long k)
{
    cout << "leftIndex NOT DEFINED!!"<<endl;
    return -1;
}
long  Operator::rightIndex(unsigned long k)
{
    cout << "leftIndex NOT DEFINED!!"<<endl;
    return -1;
}


void Operator::print()
{
    for (int i = 0; i < 6; ++i) {
        for(int j = 0; j < 12; ++j)
            cout<< setw(10) << getEntry(i,j) << " ";
     
        cout <<endl;
    }
}

Operator *Operator::operator+(Operator *a)
{
    PlusOperator *pl = new PlusOperator(this);
    pl->push_back(a);
    
    return pl;
}

Operator *Operator::operator*(Operator *a)
{
    return new TimesOperator(this,a);
}

vector<double> *Operator::operator*(vector<double> *a)
{
    vector<double> *ret = new vector<double>;
    
    unsigned long n = a->size();
    
    for (unsigned long i = 0; leftIndex(i) < n; ++i) {
        double en = 0;
        
        for (unsigned long j = leftIndex(i); j < n && j <= rightIndex(i); ++j) {
            en += getEntry(i,j)*(*a)[j];
        }
        
        
        ret->push_back(en);
    }
    
    return ret;
}


Operator *Operator::operator*(double c)
{
    return new ConstantTimesOperator(c,this);
}


PlusOperator::PlusOperator(Operator *adder)
{
    summands = new vector<Operator *>;
    summands->push_back(adder);
}

PlusOperator::~PlusOperator()
{
    //TODO: Delete entries?
    delete summands;
}

void PlusOperator::push_back(Operator *add)
{
    summands->push_back(add);
}

double PlusOperator::getEntry(unsigned long k, unsigned long j)
{
    double ret = 0;
    //    cout << "createRow " << k <<": \n";
    for(Operator *i : *summands) {
        ret += i->getEntry(k, j);
    }
    
    //    ret->print();
    
    //    cout << "end createRow" <<endl;
    
    return ret;
}

long PlusOperator::leftIndex(unsigned long row)
{
    //TODO: change -1
    long ret = -1000;
    
    for(Operator *i : *summands) {
        long li = i->leftIndex(row);
        
        if (ret == -1000)
            ret = li;
        else
            ret = min(ret, li);
    }
    
    return ret;
}

long PlusOperator::rightIndex(unsigned long row)
{
    //TODO: change -1
    long ret = -1000;
    
    for(Operator *i : *summands) {
        long li = i->rightIndex(row);
        
        ret = max(ret, li);
    }
    
    return ret;
}


TimesOperator::TimesOperator(Operator *aa, Operator *bb)
{
    a = aa;
    b = bb;
}




long TimesOperator::leftIndex(unsigned long row)
{
    
    return b->leftIndex(a->leftIndex(row));
}

long TimesOperator::rightIndex(unsigned long row)
{
    return b->rightIndex(a->rightIndex(row));
}
double TimesOperator::getEntry(unsigned long k, unsigned long j)
{
    double ret = 0;

    for (long r = a->leftIndex(k); r <= a->rightIndex(k); ++r)
        ret += a->getEntry(k,r)*b->getEntry(r,j);
    
    return ret;
}


ConstantTimesOperator::ConstantTimesOperator(double aa, Operator *bb)
{
    a = aa;
    b = bb;
}




long ConstantTimesOperator::leftIndex(unsigned long row)
{
    return b->leftIndex(row);
}

long ConstantTimesOperator::rightIndex(unsigned long row)
{
    return b->rightIndex(row);
}
double ConstantTimesOperator::getEntry(unsigned long k, unsigned long j)
{
    return a*b->getEntry(k, j);
}



ShiftOperator::ShiftOperator(Operator *a, long s)
{
    adder = a;
    shift = s;
}

double ShiftOperator::getEntry(unsigned long row, unsigned long col)
{
    return adder->getEntry(row + shift, col);
}


long ShiftOperator::leftIndex(unsigned long row)
{
    return adder->leftIndex(row + shift);
}

long ShiftOperator::rightIndex(unsigned long row)
{
    return adder->rightIndex(row + shift);
}


SkipOperator::SkipOperator(Operator *a, unsigned long rst,unsigned long rskip,unsigned long cst,unsigned long cskip)
{
    op=a;
    row_start = rst;
    row_skip = rskip;
    col_start = cst;
    col_skip = cskip;
}

unsigned long SkipOperator::shiftRow(unsigned long r)
{
    return row_skip*r + row_start;
}

unsigned long SkipOperator::shiftColumn(unsigned long r)
{
    return col_skip*r + col_start;
    
}



double SkipOperator::getEntry(unsigned long row, unsigned long col)
{
    return op->getEntry(shiftRow(row), shiftColumn(col));
}


long SkipOperator::leftIndex(unsigned long row)
{
    return (op->leftIndex(shiftRow(row))-col_start)/col_skip;
}

long SkipOperator::rightIndex(unsigned long row)
{
    return (op->rightIndex(shiftRow(row))-col_start)/col_skip;
}







