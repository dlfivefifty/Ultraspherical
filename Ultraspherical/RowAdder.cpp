//
//  RowAdder.c
//  Ultraspherical
//
//  Created by Sheehan Olver on 8/10/13.
//  Copyright (c) 2013 School of Mathematics and Statistics, The University of Sydney. All rights reserved.
//

#include <stdio.h>

#include "RowAdder.h"

RowAdder::RowAdder()
{
    
}

double RowAdder::getEntry(long k, long j)
{
    cout << "getEntry NOT DEFINED!!"<<endl;
    return NULL;
}

long RowAdder::leftIndex(unsigned long k)
{
    cout << "leftIndex NOT DEFINED!!"<<endl;
    return NULL;
}
long  RowAdder::rightIndex(unsigned long k)
{
    cout << "leftIndex NOT DEFINED!!"<<endl;
    return NULL;
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


TimesRowAdder::TimesRowAdder(RowAdder *adder)
{
    summands = new vector<RowAdder *>;
    summands->push_back(adder);
}

void TimesRowAdder::push_back(RowAdder *add)
{
    summands->push_back(add);
}



long TimesRowAdder::leftIndex(unsigned long row)
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

long TimesRowAdder::rightIndex(unsigned long row)
{
    //TODO: change -1
    long ret = -1000;
    
    for(RowAdder *i : *summands) {
        long li = i->rightIndex(row);
        
        ret = max(ret, li);
    }
    
    return ret;
}
double TimesRowAdder::getEntry(long k, long j)
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




double DerivativeRowAdder::getEntry(long row, long col)
{
    if (row == col-2)
        return (4 + 2*(row));
    else
        return 0;
}


long DerivativeRowAdder::leftIndex(unsigned long row)
{
    return row+2;
}

long DerivativeRowAdder::rightIndex(unsigned long row)
{
    return row+2;
}




double ConversionRowAdder::getEntry(long row, long col)
{
    double c = (double) col;
    
    if (  col == 0 && row ==0)
        return 1;
    else if (row == col)
        return 1/(2.*(c+1));
    else if (row +2 == col)
        return -1/(2*(c+1))-1/(2*(c-1));
    else if (row + 4 == col)
        return 1/(2*(c-1));
    else
        return 0;
}


long ConversionRowAdder::leftIndex(unsigned long row)
{
    return row;
}

long ConversionRowAdder::rightIndex(unsigned long row)
{
    return row + 4;
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
