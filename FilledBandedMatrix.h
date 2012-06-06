//
//  BandedMatrix.h
//  Ultraspherical
//
//  Created by Sheehan Olver on 04/06/2012.
//  Copyright (c) 2012 School of Mathematics and Statistics, The University of Sydney. All rights reserved.
//

#ifndef Ultraspherical_BandedMatrix_h
#define Ultraspherical_BandedMatrix_h

#include <iostream>
#include <vector>

using namespace std;

class FilledRow
{
    vector<double> entries;
    double fill;  
    int index;
public:
    FilledRow(int);
    int size();
    double operator[] (int);
    void setEntry(int,double,bool);      
    void setEntry(int,double);  
    void push_back(double);
    void setFill(double);
    void increaseSize();       
    void increaseSize(int);           
    double getFill();
    int rightIndex();
    int leftIndex();
};




class FilledBandedMatrix
{
    vector<FilledRow> rows;   
    int lowerIndex;
//    int upper;
    
public:
    FilledBandedMatrix(int,int);
    int lower();
//    int upper();
//    void increaseUpper();
//    void increaseUpper(int);    
    int size();
    int columnSize();
    int columnSize(int);
    FilledRow createRow(int);
    void increaseSize();    
//    double getEntry(int,int);
    FilledRow back();    
    
    
    void setEntry(int,int,double);    
    FilledRow operator[] (int);
    void print();
    

    
    int leftIndex(int);
    int rightIndex(int);    
    
    
// apply GIvens to two rows, and vector rhs at same time    
    void applyGivens(int,int,vector<double> *);    
    
    
    double rowDot(int,int,int,vector<double> *);
};



#endif
