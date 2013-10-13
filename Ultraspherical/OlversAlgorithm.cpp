//
//  OlversAlgorithm.cpp
//  Ultraspherical
//
//  Created by Sheehan Olver on 10/10/13.
//  Copyright (c) 2013 School of Mathematics and Statistics, The University of Sydney. All rights reserved.
//

#include "OlversAlgorithm.h"




Operator *ConstantOperator(double c)
{
    double *one = (double *)malloc(sizeof(double));
    
    *one = c;
    
    return new ToeplitzOperator(one, 1,0);
}




RSOperator::RSOperator()
{
}

double RSOperator::getEntry(long row, long col)
{
    double c = (double) col + 2;
    double r = 1/((double)4 + 2*row);
    
    if (row == col + 2)
        return r/(2.*(c+1));
    else if (row +2 == col + 2)
        return r*(-1/(2*(c+1))-1/(2*(c-1)));
    else if (row + 4 == col + 2)
        return r/(2*(c-1));
    
                

    return 0;
}


long RSOperator::leftIndex(unsigned long row)
{
    if (row < 2) {
        return 0;
    } else
        return row - 2;
}

long RSOperator::rightIndex(unsigned long row)
{
    return row + 2;
}


Operator *RSListOperator(double *a, double *s, unsigned long n)
{
    if(n < 1)
        return NULL;
        
    
    PlusOperator *pl = new PlusOperator(new DoubleTimesOperator(a[0], new RSOperator()));
    pl->push_back(ConstantOperator(s[0]));
    
    
    
    return pl;
}



