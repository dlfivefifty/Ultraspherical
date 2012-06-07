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
    
    double error = 1;
    int col = -1;
    int row1;
    
    
    while (error > 10E-17) {
        
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
            
            
            
            B.applyGivens(row1,row2,&c);   
            
            error = max(error,fabs(c[row1]));
            error = max(error,fabs(c[row2]));            
        }
        
    }
    
    //	std:cout<<"R: "<<endl;
    //    B.print();
    
    //    printvec(c);       
    
    
    vector<double> r = c;
    
    r[col] = c[col]/B.getEntry(col, col);
    
    vector<double> s;
    for(int i = 0; i < B[col].fillSize(); i++)
        s.push_back(B[col].fillGenerate(i,col)*r[col]);
    int rbnd;  
    int csize = col+1;
    
    for(int row = csize - 1; row >= 0; row--)
    {
        rbnd = B.rightIndex(row);
        if(rbnd >= csize) {
            rbnd = csize-1;
            
            r[row] = (c[row]-B.rowDot(row,row+1,rbnd,&r))/B.getEntry(row,row);
        } else {      
            //scont is the contribution from higher fills
            double scont = 0;
            for(int i = 0; i < B[row].fillSize(); i++)  
                scont += s[i]*B[row].getFill(i);    
        
            r[row] = (c[row]-B.rowDot(row,row+1,rbnd,&r) - scont)/B.getEntry(row,row);
            
            for(int i = 0; i < B[row].fillSize(); i++)            
                s[i]+=B[rbnd].fillGenerate(i,rbnd)*r[rbnd];                            
        }
    }
    
    
    //    printvec(r);
    
    return r;
}