
#include "MultiplicationMatrices.h"
#include <math.h>
#include "AdaptiveQR.h"
// TODO: Better to construct the Hankel vector in reverse so that it can be pop_backed as the ConvertMult vector is formed.


ToeplitzRowAdder::ToeplitzRowAdder(vector<double> *ain)
{
    //Assume symmetric
    a = ain;
}

double ToeplitzRowAdder::getEntry(long row, long col)
{
    long ind = abs(row - col);
    
    if (ind >= a->size())
        return 0;
    else if (ind == 0)
        return (*a)[ind];
    else
        return .5*(*a)[ind];
}


long ToeplitzRowAdder::leftIndex(unsigned long row)
{
    if (row + 1 < a->size())
        return 0;
    else
        return row+1-a->size();
}

long ToeplitzRowAdder::rightIndex(unsigned long row)
{
    return row-1+a->size();
}


HankelRowAdder::HankelRowAdder(vector<double> *ain)
{
    //Assume symmetric
    a = ain;
}

double HankelRowAdder::getEntry(long row, long col)
{
    if (row < 0)
        return 0;
    
    long ind = row + col;
    
    if (ind >= a->size())
        return 0;
    else
        return .5*(*a)[ind];
}


long HankelRowAdder::leftIndex(unsigned long row)
{
    if (row < a->size())
        return 0;
    else
        return row;
}

long HankelRowAdder::rightIndex(unsigned long row)
{
    if (row <= a->size())
        return a->size() - 1;
    else
        return row;
    
}


RowAdder *MultiplicationRowAdder(vector<double> *args)
{
    vector<double> *ha = new vector<double>;
    
    vector<double>::iterator it = args->begin();
    
    
    for (    ++it; it < args->end(); ++it) {
        ha->push_back(*it);
    }
    
    PlusRowAdder *pl = new PlusRowAdder(new ToeplitzRowAdder(args));
    pl->push_back(new ShiftRowAdder(new HankelRowAdder(ha),-1));
    
    return pl;
}





FilledBandedMatrix *DirichletD2ConvertMultiplicationMatrix(vector<double> *a)
{
    PlusRowAdder *pl =
    new   PlusRowAdder(new DerivativeRowAdder());
    pl->push_back(new TimesRowAdder(new ConversionRowAdder(),MultiplicationRowAdder(a)));
    
    
    FilledBandedMatrix *ret = new FilledBandedMatrix(-1- (int)a->size(),new ShiftRowAdder(pl,-2));
    
    FilledRow *drrow = new FilledRow(0,RowFiller::dirichlet(1,0));
    
    drrow->increaseSize();
    ret->push_back(drrow);
    
    drrow = new FilledRow(0,RowFiller::dirichlet(0,1));
    drrow->increaseSize();
    ret->push_back(drrow);
    
    return ret;
}








