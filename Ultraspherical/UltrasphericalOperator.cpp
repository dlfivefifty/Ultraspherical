
#include "UltrasphericalOperator.h"
#include <cmath>
#include "AdaptiveQR.h"
// TODO: Better to construct the Hankel vector in reverse so that it can be pop_backed as the ConvertMult vector is formed.



DerivativeOperator::DerivativeOperator(unsigned int l, unsigned int m)
{
    from = l;
    to = m;
}

long factorial(long x, long result = 1) {
    if (x == 1 || x ==0) return 1; else return factorial(x - 1, x * result);
}

double DerivativeOperator::getEntry(unsigned long row, unsigned long col)
{
    if (from == 0) {
        if (row == col - to)
            return pow((double)2,(double)(to-1))*factorial(to-1)*(to + row);
        else
            return 0;

    } else {
        if (row == col + from - to)
            return (4 + 2*(row));
        else
            return 0;
    }
}


long DerivativeOperator::leftIndex(unsigned long row)
{
    return row+2;
}

long DerivativeOperator::rightIndex(unsigned long row)
{
    return row+2;
}


ConversionOperator::ConversionOperator(unsigned int l, unsigned int m)
{
    from = l;
    to = m;
}

double ConversionOperator::getEntry(unsigned long row, unsigned long col)
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


long ConversionOperator::leftIndex(unsigned long row)
{
    return row;
}

long ConversionOperator::rightIndex(unsigned long row)
{
    return row + (to - from)*2;
}


ToeplitzOperator::ToeplitzOperator(double *ain, unsigned long l, long ind)
{
    //Assume symmetric
    length =l;
    a = ain;
    index = ind;
}

ToeplitzOperator::ToeplitzOperator(vector<double> *ain, long ind)
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

double ToeplitzOperator::getEntry(unsigned long row, unsigned long col)
{
    long ind = row - col + index;
    
    if (ind >= length)
        return 0;
    else
        return a[ind];
}


long ToeplitzOperator::leftIndex(unsigned long row)
{
    if (row < index)
        return 0;
    else
        return row-index;
}

long ToeplitzOperator::rightIndex(unsigned long row)
{
    return row+length+index;
}


HankelOperator::HankelOperator(vector<double> *ain)
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

HankelOperator::HankelOperator(double *ain, unsigned long l)
{
    a = ain;
    length = l;
}

double HankelOperator::getEntry(unsigned long row, unsigned long col)
{    
    long ind = row + col;
    
    if (ind >= length)
        return 0;
    else
        return a[ind];
}


long HankelOperator::leftIndex(unsigned long row)
{
    if (row < length)
        return 0;
    else
        return row;
}

long HankelOperator::rightIndex(unsigned long row)
{
    if (row <= length)
        return length - 1;
    else
        return row;
    
}


Operator *MultiplicationToeplitzOperator(unsigned int lambda, vector<double> *ain)
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
        
        return new ToeplitzOperator(a, length, index);
    }
    
    throw "MultiplicationToeplitzOperator not defined";
}

Operator *MultiplicationHankelOperator(unsigned int lambda, vector<double> *ain)
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
        
        return new ShiftOperator(new HankelOperator(a, length),-1);
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
        
        return new HankelOperator(a, length);
        
    }
    
    throw "MultiplicationHankelOperator not defined";
}


Operator *MultiplicationOperator(unsigned int lambda, vector<double> *args)
{
    Operator *toep = MultiplicationToeplitzOperator(lambda,args);
    Operator *han = MultiplicationHankelOperator(lambda,args);
    
    if (han == NULL) {
        return toep;
    }
    
    
    PlusOperator *pl = new PlusOperator(toep);
    pl->push_back(han);
    
    return pl;
}





FilledBandedOperator *DirichletD2ConvertMultiplicationMatrix(vector<double> *a, vector<double> *b)
{
    PlusOperator *pl =
    new   PlusOperator(new DerivativeOperator(0,2));
    if (a != NULL)
        pl->push_back(new TimesOperator(new ConversionOperator(1,2),new TimesOperator(MultiplicationOperator(1,a),new DerivativeOperator(0,1))));
    
    if (b!= NULL)
        pl->push_back(new TimesOperator(new ConversionOperator(0,2),MultiplicationOperator(0,b)));
    
    
    FilledBandedOperator *ret = new FilledBandedOperator(-1- (int)b->size(),new ShiftOperator(pl,-2));
    
    FilledRow *drrow = new FilledRow(0,RowFiller::dirichlet(1,0));
    
    drrow->increaseSize();
    ret->push_back(drrow);
    
    drrow = new FilledRow(0,RowFiller::dirichlet(0,1));
    drrow->increaseSize();
    ret->push_back(drrow);
    
    return ret;
}








