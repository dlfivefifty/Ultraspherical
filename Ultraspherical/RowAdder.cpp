//
//  RowAdder.c
//  Ultraspherical
//
//  Created by Sheehan Olver on 8/10/13.
//  Copyright (c) 2013 School of Mathematics and Statistics, The University of Sydney. All rights reserved.
//

#include <stdio.h>

#include "RowAdder.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>



RowAdder::RowAdder()
{
    
}

double RowAdder::getEntry(long k, long j)
{
    cout << "getEntry NOT DEFINED!!"<<endl;
    return -1;
}

long RowAdder::leftIndex(unsigned long k)
{
    cout << "leftIndex NOT DEFINED!!"<<endl;
    return -1;
}
long  RowAdder::rightIndex(unsigned long k)
{
    cout << "leftIndex NOT DEFINED!!"<<endl;
    return -1;
}


void RowAdder::print()
{
    for (int i = 0; i < 6; ++i) {
        for(int j = 0; j < 12; ++j)
            cout<< setw(10) << getEntry(i,j) << " ";
     
        cout <<endl;
    }
}


PlusRowAdder::PlusRowAdder(RowAdder *adder)
{
    summands = new vector<RowAdder *>;
    summands->push_back(adder);
}

void PlusRowAdder::push_back(RowAdder *add)
{
    summands->push_back(add);
}

double PlusRowAdder::getEntry(long k, long j)
{
    double ret = 0;
    //    cout << "createRow " << k <<": \n";
    for(RowAdder *i : *summands) {
        ret += i->getEntry(k, j);
    }
    
    //    ret->print();
    
    //    cout << "end createRow" <<endl;
    
    return ret;
}

long PlusRowAdder::leftIndex(unsigned long row)
{
    //TODO: change -1
    long ret = -1000;
    
    for(RowAdder *i : *summands) {
        long li = i->leftIndex(row);
        
        if (ret == -1000)
            ret = li;
        else
            ret = min(ret, li);
    }
    
    return ret;
}

long PlusRowAdder::rightIndex(unsigned long row)
{
    //TODO: change -1
    long ret = -1000;
    
    for(RowAdder *i : *summands) {
        long li = i->rightIndex(row);
        
        ret = max(ret, li);
    }
    
    return ret;
}


TimesRowAdder::TimesRowAdder(RowAdder *aa, RowAdder *bb)
{
    a = aa;
    b = bb;
}




long TimesRowAdder::leftIndex(unsigned long row)
{
    
    return b->leftIndex(a->leftIndex(row));
}

long TimesRowAdder::rightIndex(unsigned long row)
{
    return b->rightIndex(a->rightIndex(row));
}
double TimesRowAdder::getEntry(long k, long j)
{
    double ret = 0;

    for (long r = a->leftIndex(k); r <= a->rightIndex(k); ++r)
        ret += a->getEntry(k,r)*b->getEntry(r,j);
    
    return ret;
}


DoubleTimesRowAdder::DoubleTimesRowAdder(double aa, RowAdder *bb)
{
    a = aa;
    b = bb;
}




long DoubleTimesRowAdder::leftIndex(unsigned long row)
{
    return b->leftIndex(row);
}

long DoubleTimesRowAdder::rightIndex(unsigned long row)
{
    return b->rightIndex(row);
}
double DoubleTimesRowAdder::getEntry(long k, long j)
{
    return a*b->getEntry(k, j);
}



ShiftRowAdder::ShiftRowAdder(RowAdder *a, long s)
{
    adder = a;
    shift = s;
}

double ShiftRowAdder::getEntry(long row, long col)
{
    return adder->getEntry(row + shift, col);
}


long ShiftRowAdder::leftIndex(unsigned long row)
{
    return adder->leftIndex(row + shift);
}

long ShiftRowAdder::rightIndex(unsigned long row)
{
    return adder->rightIndex(row + shift);
}




