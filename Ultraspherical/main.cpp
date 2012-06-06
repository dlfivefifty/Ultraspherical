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


int main(int argc, const char * argv[])
{    
    FilledBandedMatrix bnd(-1,1);
    clock_t t1; clock_t t2;
    t1 = clock();
    
    bnd.increaseSize();
    bnd.increaseSize();
    bnd.increaseSize();
    bnd.increaseSize();
    bnd.increaseSize();
    bnd.increaseSize();
    
    bnd.increaseUpper();
    
//    cout<<"B: "<<endl;
//        bnd.print();
//    
//    cout<<"GB: "<<endl;    
//    bnd.applyGivens(0,1,0,1);
//    bnd.print();
    
    
//    
//    
    vector<double> b;
    
    for(long i = 0; i < 200000; i++)
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

