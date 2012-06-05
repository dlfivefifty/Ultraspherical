//
//  AdaptiveQR.cpp
//  Ultraspherical
//
//  Created by Sheehan Olver on 5/06/12.
//  Copyright (c) 2012 School of Mathematics and Statistics, The University of Sydney. All rights reserved.
//

#include "AdaptiveQR.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


void printvec(vector<double> c)
{
    cout <<"c: ";
    vector<double>::iterator it;
    
    for(it = c.begin(); it < c.end(); it++)
    {
        cout << *it << ", ";
    }
    
    cout << endl;
}



vector<double> QRSolve(FilledBandedMatrix B,vector<double> c)
{  
    B.increaseUpper(-B.lower());
    
    double error = 1;
    int col = -1;
    int row1;
    
    printvec(c);
    
    while (error > 10E-17) {
        cout<<"B: "<<endl;B.print();        
    printvec(c);   

        
        col++;
        row1 = col;
        
        error = fabs(c[row1]);
        
        for(int row2 = row1 + 1; row2 < B.columnSize(col); row2++)
        {
            
            if(row2 >= c.size()) {
                c.push_back(0);
            }
            if(row2 >= B.size()) {
                B.increaseSize();
            }
            
            double a = B.getEntry(row1, col);
            double b = B.getEntry(row2, col);
            double norm = sqrt(a*a + b*b);
            
            a = a/norm;
            b = b/norm;
            
            double c1 = c[row1];
            double c2 = c[row2];
            
            c[row1] = a*c1 + b*c2;
            c[row2] = -b*c1 + a*c2;
            
            error = max(error,fabs(c[row1]));
            error = max(error,fabs(c[row2]));
            
            
            B.applyGivens(row1,row2,a,b);   
        }

    }
    
	std:cout<<"R: "<<endl;
    B.print();
    
    printvec(c);       
    
    
    vector<double> r = c;
    
    r[col] = c[col]/B.getEntry(col, col);
    double s = r[col];
    int rbnd;    
    for(int row = col - 1; row >= 0; row--)
    {
        rbnd = B.rightIndex(row);
        if(rbnd >= c.size()) {
            rbnd = c.size()-1;
            
            r[row] = (c[row]-B.rowDot(row,row+1,rbnd,r))/B.getEntry(row,row);
        } else {            
            r[row] = (c[row]-B.rowDot(row,row+1,rbnd,r) - s*B.getEntry(row,rbnd+1))/B.getEntry(row,row);
            
            s+=r[rbnd];
        }
    }
    
    printvec(r);
    
    return c;
}