//
//  AdaptiveQR.cpp
//  Ultraspherical
//
//  Created by Sheehan Olver on 5/06/12.
//  Copyright (c) 2012 School of Mathematics and Statistics, The University of Sydney. All rights reserved.
//

#include "AdaptiveQR.h"

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <iostream>


void printvec(vector<double> c)
{
//    cout <<"c: ";

    
    for(long i = 0; i < c.size()-1; i++)
    {
        cout << setprecision(16)<< c[i] << ", ";
    }
    
    cout << setprecision(16)<< c[c.size()-1] << endl;    
}



vector<double> *QRSolve(FilledBandedOperator *B,vector<double> c)
{  
    
    double error = 1;
    long col = -1;
    long row1;
    long cn = c.size();
    
    
//    cout<<"colsize"<<B->columnSize(0)<<endl;
    
    //TODO:better error
    while (error > 1E-100  || row1 < cn) {
        
        col++;
        row1 = col;
        
        error = 0;
        
        //fabs(c[row1]);
        
        //    cout<< "col: " << col << " colsize: "<<B->columnSize(col)<<endl;        
        
//printvec(c);
        
        for(long row2 = row1 + 1; row2 < B->columnSize(col); row2++)
        {

            
            if(row2 >= c.size()) {
                c.push_back(0);
            }
            if(row2 >= B->size()) {
                B->increaseSize();
            }
            
            
            
            B->applyGivens(row1,row2,&c);   
            
            error = max(error,fabs(c[row1]/B->getEntry(row1,row1)));
            error = max(error,fabs(c[row2]/B->getEntry(row2,row2)));
        }
        
    }
    

    
    //Back substitition:
    
    
    vector<double> *r = new vector<double>;
    *r = c;
    
    
//        printvec(c);    
    
    (*r)[col] = c[col]/B->getEntry(col, col,true);
    
    vector<double> s;
    for(int i = 0; i < (*B)[col]->fillSize(); i++)
        s.push_back((*B)[col]->fillGenerate(i,col)*(*r)[col]);
    unsigned long rbnd;
    long csize = col+1;
    
    for(long row = csize - 1; row >= 0; row--)
    {
        rbnd = B->rightIndex(row);
        if(rbnd >= csize) {
            rbnd = csize-1;
            
            (*r)[row] = (c[row]-B->rowDot(row,row+1,rbnd,r))/B->getEntry(row,row,true);
        } else {      
            //scont is the contribution from higher fills
            double scont = 0;
            for(int i = 0; i < (*B)[row]->fillSize(); i++)
                scont += s[i]*(*B)[row]->getFill(i);    
        
            (*r)[row] = (c[row]-B->rowDot(row,row+1,rbnd,r) - scont)/B->getEntry(row,row,true);
//            
            for(int i = 0; i < (*B)[row]->fillSize(); i++)
            {
//                (*B)[rbnd];
                s[i]+=(*B)[rbnd]->fillGenerate(i,rbnd)*(*r)[rbnd];
            }
        }
    }
    
    
    //    printvec(r);
    
    return r;
}
