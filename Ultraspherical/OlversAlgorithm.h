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
#include "UltrasphericalOperator.h"
#include "FilledBandedOperator.h"

Operator *ConstantOperator(double);

class RSOperator : public Operator
{
    
public:
    RSOperator();
    
    virtual double getEntry(unsigned long,unsigned long);
    virtual long leftIndex(unsigned long);
    virtual long rightIndex(unsigned long);
};

Operator *RSListOperator(double *, double *, unsigned long n);


vector<vector<double> *> *adaptiveSylvester(Operator *Sin, vector<double> *f);
vector<vector<double> *> *poisson(vector<double> *);


#endif /* defined(__Ultraspherical__OlversAlgorithm__) */
