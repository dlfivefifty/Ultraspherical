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
    void setEntry(int,double);  
    void setEntry(int,double,bool);      
    void push_back(double);
    void setFill(double);
    void increaseSize();       
    void increaseSize(int);     
    void dropFirst();  
    double getFill();
    int leftIndex();
    int rightIndex();    
};




class FilledBandedMatrix
{
    vector<FilledRow> rows;
    int lowerIndex;
    
    
public:
    FilledBandedMatrix(int);
    int lower();
    int upper();
    void dropFirst(int row);
    int size();
    int columnSize();
    int columnSize(int);
    FilledRow createRow(int);
    void increaseSize();    
    double getEntry(int,int);
    void setEntry(int,int,double);    
    void setEntry(int,int,double,bool);        
    FilledRow operator[] (int);
    void print();
    
    int leftIndex(int);
    int rightIndex(int);    
    
    
    void applyGivens(int,int,vector<double> *);    
    
    
    double rowDot(int,int,int,vector<double> *);
};



#endif
