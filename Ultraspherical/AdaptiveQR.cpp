//
//  AdaptiveQR.cpp
//  Ultraspherical
//
//  Created by Sheehan Olver on 5/06/12.
//  Copyright (c) 2012 School of Mathematics and Statistics, The University of Sydney. All rights reserved.
//

#include "AdaptiveQR.h"

#include <math.h>




vector<double> QRSolve(FilledBandedMatrix B,vector<double> c)
{  
    B.increaseUpper(-B.lower());
    
    double error = 1;
    int col = -1;
    int row1;
    
    
    while (error > 10E-17) {
        col++;
        row1 = col;
        
        error = 0;
        
        for(int row2 = row1 + 1; row2 < B.columnSize(col); row2++)
        {
            
            if(row2 >= c.size()) {
                c.push_back(0);
            }
            
            
            double a = B.getEntry(row1, col);
            double b = B.getEntry(row2, col);
            double norm = sqrt(a*a + b*b);
            
            a = a/norm;
            b = b/norm;
            
            
            c[row1] = a*c[row1] + b*c[row2];
            c[row2] = -b*c[row1] + a*c[row2];
            
            error = max(error,c[row1]);
            error = max(error,c[row2]);
            
            
            B.applyGivens(row1,row2,a,b);            
        }

    }
    
std:cout<<"R: "<<endl;
    B.print();
    
    return c;
}