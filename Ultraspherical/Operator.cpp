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
#include <math.h>



Operator::Operator()
{
    
}

double Operator::getEntry(long k, long j)
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


PlusOperator::PlusOperator(Operator *adder)
{
    summands = new vector<Operator *>;
    summands->push_back(adder);
}

void PlusOperator::push_back(Operator *add)
{
    summands->push_back(add);
}

double PlusOperator::getEntry(long k, long j)
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
double TimesOperator::getEntry(long k, long j)
{
    double ret = 0;

    for (long r = a->leftIndex(k); r <= a->rightIndex(k); ++r)
        ret += a->getEntry(k,r)*b->getEntry(r,j);
    
    return ret;
}


DoubleTimesOperator::DoubleTimesOperator(double aa, Operator *bb)
{
    a = aa;
    b = bb;
}




long DoubleTimesOperator::leftIndex(unsigned long row)
{
    return b->leftIndex(row);
}

long DoubleTimesOperator::rightIndex(unsigned long row)
{
    return b->rightIndex(row);
}
double DoubleTimesOperator::getEntry(long k, long j)
{
    return a*b->getEntry(k, j);
}



ShiftOperator::ShiftOperator(Operator *a, long s)
{
    adder = a;
    shift = s;
}

double ShiftOperator::getEntry(long row, long col)
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




