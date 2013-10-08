
#include "UltrasphericalRowAdder.h"
#include <math.h>
#include "AdaptiveQR.h"
// TODO: Better to construct the Hankel vector in reverse so that it can be pop_backed as the ConvertMult vector is formed.



DerivativeRowAdder::DerivativeRowAdder(unsigned int l, unsigned int m)
{
    from = l;
    to = m;
}

long factorial(long x, long result = 1) {
    if (x == 1) return result; else return factorial(x - 1, x * result);
}

double DerivativeRowAdder::getEntry(long row, long col)
{
    if (from == 0) {
        if (row == col - to)
            return pow(2,to-1)*factorial(to-1)*(to + row);
        else
            return 0;

    } else {
        if (row == col + from - to)
            return (4 + 2*(row));
        else
            return 0;
    }
}


long DerivativeRowAdder::leftIndex(unsigned long row)
{
    return row+2;
}

long DerivativeRowAdder::rightIndex(unsigned long row)
{
    return row+2;
}


ConversionRowAdder::ConversionRowAdder(unsigned int l, unsigned int m)
{
    from = l;
    to = m;
}

double ConversionRowAdder::getEntry(long row, long col)
{
    double c = (double) col;
    
    if (from == to && row == col) {
        return 1;
    } else if (from == 0) {
        switch (to) {
            case 1:
                if (row == 0 && col == 0)
                    return 1;
                else if (row == col)
                    return .5;
                else if (row +2 == col)
                    return -.5;
                
                break;
                
                
            case 2:
                
                if (  col == 0 && row ==0)
                    return 1;
                else if (row == col)
                    return 1/(2.*(c+1));
                else if (row +2 == col)
                    return -1/(2*(c+1))-1/(2*(c-1));
                else if (row + 4 == col)
                    return 1/(2*(c-1));
                
                
                break;
                
            default:
                break;
        }
    } else if (from == 1) {
        switch (to) {
            case 2:
                if (row == col)
                    return 1/((double) (col+1));
                else if (row +2 == col)
                    return -1/((double) (col+1));
                    
                break;
                
            default:
                break;
        }
    }


    return 0;
}


long ConversionRowAdder::leftIndex(unsigned long row)
{
    return row;
}

long ConversionRowAdder::rightIndex(unsigned long row)
{
    return row + (to - from)*2;
}



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


RowAdder *MultiplicationRowAdder(unsigned int lambda, vector<double> *args)
{
    if (lambda ==0) {
        vector<double> *ha = new vector<double>;
        
        vector<double>::iterator it = args->begin();
        
        
        for (    ++it; it < args->end(); ++it) {
            ha->push_back(*it);
        }
        
        PlusRowAdder *pl = new PlusRowAdder(new ToeplitzRowAdder(args));
        pl->push_back(new ShiftRowAdder(new HankelRowAdder(ha),-1));
        
        return pl;
    } else {
        //TODO: Implement
        throw "Higher order multiplication not implemented";
    }
}





FilledBandedMatrix *DirichletD2ConvertMultiplicationMatrix(vector<double> *a)
{
    PlusRowAdder *pl =
    new   PlusRowAdder(new DerivativeRowAdder(0,2));
    pl->push_back(new TimesRowAdder(new ConversionRowAdder(0,2),MultiplicationRowAdder(0,a)));
    
    
    FilledBandedMatrix *ret = new FilledBandedMatrix(-1- (int)a->size(),new ShiftRowAdder(pl,-2));
    
    FilledRow *drrow = new FilledRow(0,RowFiller::dirichlet(1,0));
    
    drrow->increaseSize();
    ret->push_back(drrow);
    
    drrow = new FilledRow(0,RowFiller::dirichlet(0,1));
    drrow->increaseSize();
    ret->push_back(drrow);
    
    return ret;
}








