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


//FilledRow defaultRowCreator(int k)
//{
//    
//    double w = 1;
//    
//    vector<RowFiller> fillers;
//    fillers.push_back(RowFiller::leftDirichlet());
//    fillers.push_back(RowFiller::rightDirichlet());    
//    
//    switch(k)
//    {
//        case 0:
//        {
//            FilledRow newrow(0,RowFiller::dirichlet(1,0));
//            newrow.increaseSize();//must have an element in first row to use Givens                        
//            return newrow;
//        }
//        case 1:
//        {
//            FilledRow newrow(0,RowFiller::dirichlet(0,1));
//            newrow.increaseSize();//must have an element in first row to delete            
//            return newrow;
//        }            
//        case 2:
//        {
//            FilledRow newrow(0,RowFiller::dirichlet(0,0),3);               
//            newrow.setEntry(0,w);
//            newrow.setEntry(1,1);
//            newrow.setEntry(2,-.5*w);
//        
//            return newrow;
//        }
//        default:
//        {
//            FilledRow newrow(k-2,RowFiller::dirichlet(0,0),3);       
//            
//            newrow.setEntry(k-2,.5*w);
//            newrow.setEntry(k-1,k-1);
//            newrow.setEntry(k,-.5*w);            
//            
//            return newrow;
//        }
//    }    
//}


FilledRow D2Dirichlet(int k)
{
    
//    double w = 1;
    
    switch(k)
    {
        case 0://LeftDir
        {
            FilledRow newrow(0,RowFiller::dirichlet(1,0));
            newrow.increaseSize();//must have an element in first row to use Givens                        
            return newrow;
        }
        case 1://RightDir
        {
            FilledRow newrow(0,RowFiller::dirichlet(0,1));
            newrow.increaseSize();//must have an element in first row to delete            
            return newrow;
        }     
        case 2:
        {
            FilledRow newrow(0,RowFiller::dirichlet(0,0));
            newrow.push_back(1);
            newrow.push_back(0);            
            newrow.push_back(4.-2./3.);            
            newrow.push_back(0);            
            newrow.push_back(1./6);                        
            return newrow;    
        }
        default:
        {
            FilledRow newrow(k-2,RowFiller::dirichlet(0,0));
            newrow.push_back(1./(2.*(k-1.)));
            newrow.push_back(0);            
            newrow.push_back(2.*k-k/(2.*(k-1)+(k-1)*(k-1)));            
            newrow.push_back(0);            
            newrow.push_back(1./(2.*(k+1)));                   
            return newrow;  
        }
    }    
}


int main(int argc, const char * argv[])
{    
    //First entry is lower row, it's important!
    FilledBandedMatrix bnd(-2,&D2Dirichlet);
    clock_t t1; clock_t t2;

    
    
    
    
    
        bnd.increaseSize();
        bnd.increaseSize();
        bnd.increaseSize();    
        bnd.increaseSize();    
        bnd.increaseSize(); 
        bnd.increaseSize();     
    
    
    cout<<"B: "<<endl;
    bnd.print();
    //    
//        cout<<"GB: "<<endl;    
//        bnd.applyGivens(0,1,0,1);
//        bnd.print();
//    
    
    //    
    //    
    vector<double> b;

    b.push_back(1);
    
   vector<double> c = QRSolve(bnd,b);    
    
    printvec(c);
//    
    for(long i = 0; i < 50000; i++)
        b.push_back(1);
    
    t1 = clock();
    
    c = QRSolve(bnd,b);
    
    t2 = clock();
    float tottime = ((float)(t2-t1)/CLOCKS_PER_SEC);
    cout<<"time: "<<tottime<<endl;    
    cout<<"size: "<<c.size()<<endl;
    
    cout<<"time/size: "<<tottime*300000/b.size()<<endl;    
    
    
    return 0;
}

