//
//  main.cpp
//  Ultraspherical
//
//  Created by Sheehan Olver on 04/06/2012.
//  Copyright (c) 2012 School of Mathematics and Statistics, The University of Sydney. All rights reserved.
//


#include "FilledBandedMatrix.h"

#include "AdaptiveQR.h"



int main(int argc, const char * argv[])
{    
    FilledBandedMatrix bnd(-1,1);
    
    bnd.increaseSize();
    bnd.increaseSize();
    bnd.increaseSize();
    bnd.increaseSize();
    bnd.increaseSize();
    bnd.increaseSize();
    
        bnd.print();
    
    vector<double> b;
    b.push_back(1);
    vector<double> c = QRSolve(bnd,b);
    
    
    cout<<"b "<<b.size()<<endl;
    cout<<"c "<<c.size()<<endl;
    
    
    cout<<bnd.size()<<endl;
    
    
    
    return 0;
}

