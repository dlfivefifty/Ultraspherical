//
//  OlversAlgorithm.h
//  Ultraspherical
//
//  Created by Sheehan Olver on 10/10/13.
//  Copyright (c) 2013 School of Mathematics and Statistics, The University of Sydney. All rights reserved.
//

#ifndef __Ultraspherical__OlversAlgorithm__
#define __Ultraspherical__OlversAlgorithm__

#include <iostream>
#include "UltrasphericalRowAdder.h"

RowAdder *ConstantOperator(double);

class RSRowAdder : public RowAdder
{
    
public:
    RSRowAdder();
    
    virtual double getEntry(long,long);
    virtual long leftIndex(unsigned long);
    virtual long rightIndex(unsigned long);
};

RowAdder *RSListRowAdder(double *, double *, unsigned long n);


#endif /* defined(__Ultraspherical__OlversAlgorithm__) */
