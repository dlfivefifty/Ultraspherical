
#include "MultiplicationMatrices.h"
#include <math.h>
#include "AdaptiveQR.h"
// TODO: Better to construct the Hankel vector in reverse so that it can be pop_backed as the ConvertMult vector is formed.





FilledBandedMatrix *DirichletD2ConvertMultiplicationMatrix(vector<double> a) 
{
    PlusRowAdder *pl = new PlusRowAdder(new ConversionRowAdder());
    pl->push_back(new DerivativeRowAdder());
    FilledBandedMatrix *ret = new FilledBandedMatrix(-2,new ShiftRowAdder(pl,-2));
    
    FilledRow *drrow = new FilledRow(0,RowFiller::dirichlet(1,0));
    
    drrow->increaseSize();
    ret->push_back(drrow);
    
    drrow = new FilledRow(0,RowFiller::dirichlet(0,1));
    drrow->increaseSize();
    ret->push_back(drrow);
    
    return ret;
}








