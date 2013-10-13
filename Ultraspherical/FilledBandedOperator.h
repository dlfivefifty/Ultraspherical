//
//  BandedOperator.h
//  Ultraspherical
//
//  Created by Sheehan Olver on 04/06/2012.
//  Copyright (c) 2012 School of Mathematics and Statistics, The University of Sydney. All rights reserved.
//

#ifndef Ultraspherical_BandedOperator_extras_h
#define Ultraspherical_BandedOperator_extras_h

#include <iostream>
#include <vector>

#include "Operator.h"


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



class SavedOperator : public Operator
{
public:
    Operator *adder;
    vector<FilledRow *> rows;
    
    SavedOperator(Operator *);
    ~SavedOperator();
    
    FilledRow *createRow(unsigned long);
    unsigned long size();    
    void increaseSize();
    void increaseSize(unsigned long);
    double getEntry(unsigned long,unsigned long);
    double getEntry(unsigned long,unsigned long,bool);
    
    FilledRow *getRow(unsigned long,bool);
    FilledRow *operator[] (unsigned long);
    
    long leftIndex(unsigned long);
    long rightIndex(unsigned long);
    void print();
    
    virtual Operator *operator+(Operator *);
    virtual Operator *operator*(double);
};




class FilledBandedOperator : public SavedOperator
{
    int lowerIndex;
    
    
public:
    FilledBandedOperator(const int, Operator *add);
    int lower();
    void setLower(int);
    void dropFirst(unsigned long);
    unsigned long columnSize();
    unsigned long columnSize(unsigned long);
    void setEntry(unsigned long,unsigned long,double);
    void setEntry(unsigned long,unsigned long,double,bool);
    
    void push_back(FilledRow *);
    
    void setAdder(Operator *);
    Operator *getAdder();
    
    void applyGivens(unsigned long,unsigned long,vector<double> *);
    
    
    double rowDot(unsigned long,unsigned long,unsigned long,vector<double> *);
};





#endif
