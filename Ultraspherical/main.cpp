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
    
    if (k == 0) {
        FilledRow newrow(0);        
        newrow.setFill(1);
        return newrow;
    }
    else if (k == 1) {
        FilledRow newrow(0);               
        newrow.push_back(w);
        newrow.push_back(1);
        newrow.push_back(-.5*w);
        newrow.setFill(0);
        
        return newrow;
    } else {
        FilledRow newrow(k-1);       
        
        newrow.push_back(.5*w);
        newrow.push_back(k);
        newrow.push_back(-.5*w);       
        newrow.setFill(0);        
        
        return newrow;
    }    
}



int main(int argc, const char * argv[])
{    
    FilledBandedMatrix bnd(-1,&defaultRowCreator);
    clock_t t1; clock_t t2;
    t1 = clock();
    
    
    
    
    
    //    bnd.increaseSize();
    //    bnd.increaseSize();
    //    bnd.increaseSize();    
    
    
    cout<<"B: "<<endl;
    bnd.print();
    //    
    //    cout<<"GB: "<<endl;    
    //    bnd.applyGivens(0,1,0,1);
    //    bnd.print();
    
    
    //    
    //    
    vector<double> b;
    
    for(long i = 0; i < 1; i++)
        b.push_back(1);
    
    vector<double> c = QRSolve(bnd,b);
    
    t2 = clock();
    float tottime = ((float)(t2-t1)/CLOCKS_PER_SEC);
    cout<<"time: "<<tottime<<endl;    
    cout<<"size: "<<c.size()<<endl;
    
    cout<<"time/size: "<<tottime*300000/b.size()<<endl;    
    printvec(c);
    
    
    return 0;
}

