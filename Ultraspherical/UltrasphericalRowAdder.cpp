
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


ToeplitzRowAdder::ToeplitzRowAdder(double *ain, unsigned long l, long ind)
{
    //Assume symmetric
    length =l;
    a = ain;
    index = ind;
}

ToeplitzRowAdder::ToeplitzRowAdder(vector<double> *ain, long ind)
{
    //Assume symmetric
    length = ain->size();
    
    a = (double *)malloc(length*sizeof(double));
    index = ind;
    
    int k=0;
    for (double i : *ain) {
        a[k] = i;
        k++;
    }

}

double ToeplitzRowAdder::getEntry(long row, long col)
{
    long ind = row - col + index;
    
    if (ind >= length)
        return 0;
    else
        return a[ind];
}


long ToeplitzRowAdder::leftIndex(unsigned long row)
{
    if (row < index)
        return 0;
    else
        return row-index;
}

long ToeplitzRowAdder::rightIndex(unsigned long row)
{
    return row+length+index;
}


HankelRowAdder::HankelRowAdder(vector<double> *ain)
{
    //Assume symmetric
    length = ain->size();
    
    a = (double *)malloc(length*sizeof(double));
    
    int k=0;
    for (double i : *ain) {
        a[k] = i;
        k++;
    }
}

HankelRowAdder::HankelRowAdder(double *ain, unsigned long l)
{
    a = ain;
    length = l;
}

double HankelRowAdder::getEntry(long row, long col)
{
    if (row < 0)
        return 0;
    
    long ind = row + col;
    
    if (ind >= length)
        return 0;
    else
        return a[ind];
}


long HankelRowAdder::leftIndex(unsigned long row)
{
    if (row < length)
        return 0;
    else
        return row;
}

long HankelRowAdder::rightIndex(unsigned long row)
{
    if (row <= length)
        return length - 1;
    else
        return row;
    
}


RowAdder *MultiplicationToeplitzRowAdder(unsigned int lambda, vector<double> *ain)
{
    if (lambda <= 1) {
        unsigned long length = 2*ain->size() - 1;
        
        
        double *a = (double *)malloc(length*sizeof(double));
        long index = ain->size() - 1;
        
        int k=0;
        for (double i : *ain) {
            a[index - k] = a[index + k] = (k==0?1:.5)*i;
            k++;
        }
        
        return new ToeplitzRowAdder(a, length, index);
    }
    
    throw "MultiplicationToeplitzRowAdder not defined";
}

RowAdder *MultiplicationHankelRowAdder(unsigned int lambda, vector<double> *ain)
{
    if (lambda == 0) {
        if (ain->size() <= 1)
            return NULL;
        
        unsigned long length = ain->size() - 1;
        
        
        double *a = (double *)malloc(length*sizeof(double));
        
        int k=0;
        for (double i : *ain) {
            if (k > 0) {
                a[k-1] = .5*i;
            }

            k++;
        }
        
        return new ShiftRowAdder(new HankelRowAdder(a, length),-1);
    } else if (lambda == 1) {
        if (ain->size() <= 2)
            return NULL;
        
        unsigned long length = ain->size() - 2;
        
        
        double *a = (double *)malloc(length*sizeof(double));
        
        int k=0;
        for (double i : *ain) {
            if (k > 1) {
                a[k-2] = -.5*i;
            }
            
            k++;
        }
        
        return new ShiftRowAdder(new HankelRowAdder(a, length),-1);
        
    }
    
    throw "MultiplicationHankelRowAdder not defined";
}


RowAdder *MultiplicationRowAdder(unsigned int lambda, vector<double> *args)
{
    RowAdder *toep = MultiplicationToeplitzRowAdder(lambda,args);
    RowAdder *han = MultiplicationHankelRowAdder(lambda,args);
    
    if (han == NULL) {
        return toep;
    }
    
    
    PlusRowAdder *pl = new PlusRowAdder(toep);
    pl->push_back(han);
    
    return pl;
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








