//
//  main.cpp
//  Ultraspherical
//
//  Created by Sheehan Olver on 04/06/2012.
//  Copyright (c) 2012 School of Mathematics and Statistics, The University of Sydney. All rights reserved.
//


#include "FilledBandedMatrix.h"

#include "AdaptiveQR.h"
#include <time.h>


FilledRow defaultRowCreator(int k)
{
    
    double w = 1;
    
    switch(k)
    {
        case 0:
        {
            return FilledRow::leftDirichlet();
        }
        case 1:
        {
            FilledRow newrow(0,0,&alternatingFiller,3);               
            newrow.setEntry(0,w);
            newrow.setEntry(1,1);
            newrow.setEntry(2,-.5*w);
        
            return newrow;
        }
        default:
        {
            FilledRow newrow(k-1,0,&alternatingFiller,3);       
            
            newrow.setEntry(k-1,.5*w);
            newrow.setEntry(k,k);
            newrow.setEntry(k+1,-.5*w);            
            
            return newrow;
        }
    }    
}



int main(int argc, const char * argv[])
{    
    FilledBandedMatrix bnd(-1,&defaultRowCreator);
    clock_t t1; clock_t t2;
    t1 = clock();
    
    
    
    
    
        bnd.increaseSize();
        bnd.increaseSize();
        bnd.increaseSize();    
        bnd.increaseSize();    
    
    
    cout<<"B: "<<endl;
    bnd.print();
    //    
    //    cout<<"GB: "<<endl;    
    //    bnd.applyGivens(0,1,0,1);
    //    bnd.print();
    
    
    //    
    //    
    vector<double> b;
    
    for(long i = 0; i < 100000; i++)
        b.push_back(1);
    
    vector<double> c = QRSolve(bnd,b);
    
    t2 = clock();
    float tottime = ((float)(t2-t1)/CLOCKS_PER_SEC);
    cout<<"time: "<<tottime<<endl;    
    cout<<"size: "<<c.size()<<endl;
    
    cout<<"time/size: "<<tottime*300000/b.size()<<endl;    
//    printvec(c);
    
    
    return 0;
}

