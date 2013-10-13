//
//  AdaptiveQR.h
//  Ultraspherical
//
//  Created by Sheehan Olver on 5/06/12.
//  Copyright (c) 2012 School of Mathematics and Statistics, The University of Sydney. All rights reserved.
//

#ifndef Ultraspherical_AdaptiveQR_h
#define Ultraspherical_AdaptiveQR_h

#include "FilledBandedOperator.h"

vector<double> QRSolve(FilledBandedOperator *,vector<double>);

void printvec(vector<double>);


#endif

