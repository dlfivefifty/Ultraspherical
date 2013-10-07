//
//  BandedMatrix.h
//  Ultraspherical
//
//  Created by Sheehan Olver on 04/06/2012.
//  Copyright (c) 2012 School of Mathematics and Statistics, The University of Sydney. All rights reserved.
//

#ifndef Ultraspherical_BandedMatrix_extras_h
#define Ultraspherical_BandedMatrix_extras_h

#include <iostream>
#include <vector>


#define NUMFILLERS 2

using namespace std;


double oneFiller(int);
double alternatingFiller(int);



class RowFiller
{
    double fill;
    double (*fillGenerator)(unsigned long);
public:    
    RowFiller();    
    RowFiller(double, double (*)(unsigned long));
    double getFill();
    double fillGenerate(unsigned long);
    double getEntry(unsigned long);
    void setFill(double);    
    
//    static RowFiller leftDirichlet();
//    static RowFiller rightDirichlet();    
    static RowFiller** dirichlet(unsigned long,unsigned long);        
};



class FilledRow
{
    vector<double> *entries;
    RowFiller **fillers;
    long index;
    
public:
//    FilledRow(unsigned long);
//    FilledRow(unsigned long,double);
    FilledRow(long,RowFiller **);
    FilledRow(long,vector<double> *,RowFiller **);
    FilledRow(const FilledRow &);
    ~FilledRow();    
    int size();    
    double operator[] (long);
    FilledRow* operator+(FilledRow *);
    void setEntry(long,double);
    void setEntry(long,double,bool);
    void push_back(double);
    void increaseSize();       
    void increaseSize(long);
    void dropFirst();
    double getFill(int);
    void setFill(int,double);
    int fillSize();
    double fillGenerate(int,long);
    long leftIndex();
    long rightIndex();    
    void print();

    
//    static FilledRow rightDirichlet();
//    static FilledRow leftDirichlet();    
};






class FilledBandedMatrix
{
    vector<FilledRow *> rows;
    int lowerIndex;
    
    
public:
    FilledBandedMatrix(const int);
    ~FilledBandedMatrix();
    int lower();
    void setLower(int);
    void dropFirst(unsigned long);
    unsigned long size();
    unsigned long columnSize();
    unsigned long columnSize(unsigned long);
    virtual FilledRow *createRow(unsigned long);
    void increaseSize();   
    void increaseSize(unsigned long);
    double getEntry(unsigned long,unsigned long);
    double getEntry(unsigned long,unsigned long,bool);
    void setEntry(unsigned long,unsigned long,double);
    void setEntry(unsigned long,unsigned long,double,bool);
    FilledRow *getRow(unsigned long,bool);
    FilledRow *operator[] (unsigned long);
    void print();
    
    unsigned long leftIndex(unsigned long);
    unsigned long rightIndex(unsigned long);
    
    void push_back(FilledRow *);
    
    
    void applyGivens(unsigned long,unsigned long,vector<double> *);
    
    
    double rowDot(unsigned long,unsigned long,unsigned long,vector<double> *);
};





#endif
