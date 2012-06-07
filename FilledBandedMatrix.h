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


double oneFiller(int);
double alternatingFiller(int);

class RowFiller
{
    double fill;
    double (*fillGenerator)(int); 
public:    
    RowFiller();    
    RowFiller(double, double (*)(int));
    double getFill();
    double fillGenerate(int);
    double getEntry(int k);  
    void setFill(double);    
    
    static RowFiller leftDirichlet();
    static RowFiller rightDirichlet();    
    static vector<RowFiller> dirichlet(int,int);        
};


class FilledRow
{
    vector<double> entries;
    vector<RowFiller> fillers;
    int index;
    
public:
    FilledRow(int);
    FilledRow(int,double);   
    FilledRow(int,vector<RowFiller>);  
    FilledRow(int,vector<RowFiller>,int);      
    int size();
    double operator[] (int);
    void setEntry(int,double);  
    void setEntry(int,double,bool);      
    void push_back(double);
    void increaseSize();       
    void increaseSize(int);     
    void dropFirst();  
    double getFill(int);
    void setFill(int,double);    
    int fillSize();
    double fillGenerate(int,int);
    int leftIndex();
    int rightIndex();    
    
    static FilledRow rightDirichlet();
    static FilledRow leftDirichlet();    
};






class FilledBandedMatrix
{
    vector<FilledRow> rows;
    int lowerIndex;
    FilledRow (*rowGenerator)(int);
    
    
public:
    FilledBandedMatrix(int, FilledRow (*)(int));
    int lower();
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