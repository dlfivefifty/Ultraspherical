//
//  OlversAlgorithm.cpp
//  Ultraspherical
//
//  Created by Sheehan Olver on 10/10/13.
//  Copyright (c) 2013 School of Mathematics and Statistics, The University of Sydney. All rights reserved.
//

#include "OlversAlgorithm.h"




RowAdder *IdentityOperator()
{
    double *one = (double *)malloc(sizeof(double));
    
    *one = 1;
    
    return new ToeplitzRowAdder(one, 1,0);
}




RSRowAdder::RSRowAdder()
{
}

double RSRowAdder::getEntry(long row, long col)
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


long RSRowAdder::leftIndex(unsigned long row)
{
    if (row < 2) {
        return 0;
    } else
        return row - 2;
}

long RSRowAdder::rightIndex(unsigned long row)
{
    return row + 2;
}






