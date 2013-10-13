//
//  BandedOperator.cpp
//  Ultraspherical
//
//  Created by Sheehan Olver on 04/06/2012.
//  Copyright (c) 2012 School of Mathematics and Statistics, The University of Sydney. All rights reserved.
//

#include "FilledBandedOperator.h"
#include "UltrasphericalOperator.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>

#define SHIFTROW(j)(j - index)
#define FILLROWQ(j)(j - index >= size())


#define NUMFILLERS 2

double oneFiller(unsigned long k) // for Right Dirichlet
{
    return 1;
}

double alternatingFiller(unsigned long k) // for Left Dirichlet
{
    unsigned long km = (k % 2);
    return   (double)1-2*km;
}



RowFiller::RowFiller()
{
    fill = 0;
    fillGenerator = &oneFiller;
}

RowFiller::RowFiller(double fil, double (*filr)(unsigned long))   
{
    fill = fil;
    fillGenerator = filr;
}

double RowFiller::getFill()
{
    return fill;
}

void RowFiller::setFill(double val)
{
    fill = val;
}

double RowFiller::getEntry(unsigned long k)
{
   return  fill*fillGenerator(k);
}


double RowFiller::fillGenerate(unsigned long col)
{
    return fillGenerator(col);
}


//RowFiller RowFiller::leftDirichlet()
//{
//    RowFiller flr(1,&alternatingFiller);
//    return flr;
//}
//
//RowFiller RowFiller::rightDirichlet()
//{
//    RowFiller flr(1,&oneFiller);
//    return flr;
//}

RowFiller** RowFiller::dirichlet(unsigned long a, unsigned long b)
{
    RowFiller **fil = (RowFiller **) malloc(NUMFILLERS*sizeof(RowFiller*));

    fil[0] = new RowFiller(a,&alternatingFiller);
    fil[1] = new RowFiller(b,&oneFiller);    
    
    return fil;
}
//
//FilledRow FilledRow::rightDirichlet()
//{
//    RowFiller  fillers[NUMFILLERS];
//    fillers.push_back(RowFiller::rightDirichlet());
//    FilledRow newrow(0,fillers);     
//    return newrow;
//}
//
//FilledRow FilledRow::leftDirichlet()
//{
//    vector<RowFiller> fillers;
//    fillers.push_back(RowFiller::leftDirichlet());
//    FilledRow newrow(0,fillers);        
//    return newrow;
//}


//FilledRow::FilledRow(unsigned long ind)
//{
//    index = ind;
//    fillers = (RowFiller *) malloc(sizeof(RowFiller)*NUMFILLERS);    
//}


FilledRow::FilledRow(long ind,RowFiller **flr)
{
    index = ind;
    fillers = flr;
    entries = new vector<double>;    
}


FilledRow::FilledRow(long ind,vector<double> *entrs,RowFiller **flr)
{
    index = ind;
    fillers = flr;
    
    entries = entrs;
}

FilledRow::FilledRow(const FilledRow &ref)
{
    fillers = ref.fillers;    
    index = ref.index;
    entries = ref.entries;
}

FilledRow::~FilledRow()
{
    delete entries;
    for (int i = 0; i < NUMFILLERS; i++) {
        delete fillers[i];
    }
    
    free(fillers);
}

int FilledRow::size()
{
    return (int) entries->size();
}

double FilledRow::operator[](long j)
{
    //    cout<<"size(): "<<entries.size()<<endl;
    
    if (FILLROWQ(j)) {
        double ret = 0;
        for(int i = 0; i < NUMFILLERS; i++)
            ret+=fillers[i]->getEntry(j);
        
        return ret;
    }
    else if(SHIFTROW(j) < 0)
        return 0;
    else
        return (*entries)[SHIFTROW(j)];
}


void FilledRow::setEntry(long j,double val)
{
    setEntry(j,val,false);
}
void FilledRow::setEntry(long j,double val,bool increasesize)
{    
    if(SHIFTROW(j) < 0)
        throw "Filled Rows can only grow right, not left";
    
    if(increasesize == false) {
        if(FILLROWQ(j))
            throw "Set increasesize to true to change fill rows";
        
    }
    else {
        increaseSize(j);
    }
    
    (*entries)[SHIFTROW(j)] = val;
}

void FilledRow::setFill(int i,double val)
{
    fillers[i]->setFill(val);
}

double FilledRow::getFill(int i)
{
    return fillers[i]->getFill();
}

int FilledRow::fillSize()
{
    return NUMFILLERS;
}

double FilledRow::fillGenerate(int i, long col)
{
    return fillers[i]->fillGenerate(col);
}

void FilledRow::push_back(double val)
{    
    entries->push_back(val);
}

void FilledRow::increaseSize()
{
    push_back((*this)[size()]);
}

// increases Size so that row[k] is not a fill row
void FilledRow::increaseSize(long k)
{
    for(long i  = size(); i <= SHIFTROW(k); i++)
        increaseSize();
}

void FilledRow::dropFirst()
{
    if(size() == 0)
        throw "No first to drop";
    
    entries->erase(entries->begin());
    
    index++;
}


long FilledRow::leftIndex()
{
    return index;
}

long FilledRow::rightIndex()
{
    return index + size() - 1;
}


void FilledRow::print()
{
    
    for(unsigned long i = 0; i < rightIndex()+3; i++)
    {
        cout << setw(10)<<(*this)[i];
        

        
        cout << " ";
    }
    cout << endl;
}


FilledRow *FilledRow::operator+(FilledRow *row)
{
    long li = min(this->leftIndex(),row->leftIndex());
    long ri = max
    (this->rightIndex(),row->rightIndex());
    
    FilledRow *ret = new FilledRow(li,fillers);
    
    
    for (long i = li; i <= ri; ++i) {
        ret->setEntry(i, (*this)[i] + (*row)[i],true);
    }
    
    //Assuming 0 fill in
    
    return ret;
}

SavedOperator::SavedOperator(Operator *a)
{
    adder = a;
}

SavedOperator::~SavedOperator()
{
    for (unsigned long i = 0; i < rows.size(); i++) {
        delete rows[i];
    }
}



FilledBandedOperator::FilledBandedOperator(const int low, Operator *a) : SavedOperator(a)
{
    lowerIndex = low;
    
//    increaseSize();
}





int FilledBandedOperator::lower()
{
    // Returns zero if the band starts on the diagonal
    // -1 if there is one subdiagonal...
    return lowerIndex;
}

void FilledBandedOperator::setLower(int low)
{
    lowerIndex = low;
}


unsigned long SavedOperator::size()
{
    return rows.size();
}

unsigned long FilledBandedOperator::columnSize()
{
    return rows.back()->rightIndex();
}

unsigned long FilledBandedOperator::columnSize(unsigned long col)
{
    return col - lower() + 1;
}

long SavedOperator::leftIndex(unsigned long row)
{
    long li = (*this)[row]->leftIndex();
    
    if(li > 0)
        return li;
    else
        return 0;
}

long SavedOperator::rightIndex(unsigned long row)
{
    return (*this)[row]->rightIndex();
}


FilledRow *SavedOperator::getRow(unsigned long i,bool increasesize)
{    
    if(i >= size()) {
        if(increasesize) {
            increaseSize(i);
            
            return rows[i];
        } else {
            return createRow(i);
        }
    } else {
        return rows[i];
    }
}

FilledRow * SavedOperator::operator[] (unsigned long i)
{
    return getRow(i,true);
}

double SavedOperator::getEntry(unsigned long i,unsigned long j)
{    
    return (*(*this)[i])[j];
}

double SavedOperator::getEntry(unsigned long i,unsigned long j,bool increasesize)
{    
    return (*getRow(i, increasesize))[j];
}


void SavedOperator::increaseSize()
{
    rows.push_back(createRow(size()));
}

void SavedOperator::increaseSize(unsigned long i)
{
    for(unsigned long j = size(); j <= i; j++)
        increaseSize();
}


Operator *SavedOperator::operator+(Operator *a)
{
    PlusOperator *pl = new PlusOperator(this);
    pl->push_back(a);
    
    return new SavedOperator(pl);
}

Operator *SavedOperator::operator*(Operator *a)
{
    return new SavedOperator(new TimesOperator(this,a));
}

Operator *SavedOperator::operator*(double c)
{
    return new SavedOperator(new ConstantTimesOperator(c,this));
}



void FilledBandedOperator::push_back(FilledRow *row)
{
    rows.push_back(row);
}


void FilledBandedOperator::setAdder(Operator *a)
{
    adder = a;
}

Operator * FilledBandedOperator::getAdder()
{
    return adder;
}


FilledRow *SavedOperator::createRow(unsigned long k)
{
    
    
    vector<double> *newrow = new vector<double>;
    
    long li = max((long)0,adder->leftIndex(k));
    
    for (long j = li; j <= adder->rightIndex(k); ++j)
        newrow->push_back(adder->getEntry(k,j));
    
    return new FilledRow(li,newrow,RowFiller::dirichlet(0, 0));
}



void FilledBandedOperator::dropFirst(unsigned long row)
{
    rows[row]->dropFirst();
}


void FilledBandedOperator::setEntry(unsigned long i,unsigned long j,double val)
{
    setEntry(i,j,false);
}

void FilledBandedOperator::setEntry(unsigned long i,unsigned long j,double val,bool increasesize)
{    
    for(unsigned long k = size(); k <= i; k++)
    {
        rows.push_back(createRow(k));
    }
    
    rows[i]->setEntry(j,val,increasesize);    
}



void SavedOperator::print()
{

    for(unsigned long i = 0; i < size()+3; i++)
    {  
        for(unsigned long j = 0; j <= rightIndex(i); j++)
        {
            
            cout << setw(10)<<getEntry(i,j);
        }
//        cout <<"F: ";
        for(unsigned long j =  rightIndex(i) + 1; j <= rightIndex(i) + 4; j++)
        {
            cout <<setw(10) << getEntry(i,j);
        }        
        
        cout << endl;
    }
}


void FilledBandedOperator::applyGivens(unsigned long row1, unsigned long row2, vector<double> *c)
{
    unsigned long col1 = row1;    
    if(leftIndex(row2) < leftIndex(row1))
        throw "Givens should be applied left to right and top to bottom, in order";    
    if (leftIndex(row2) > col1)
        throw "Cannot change elements to the left";
    if (rightIndex(row1) > rightIndex(row2))
        throw "Probably a bug: rightIndex(row1) > rightIndex(row2)";
    if (row2 > size())
        throw "Probably a bug: row2 > size";
    
    increaseSize(max(row1,row2));
    
    
    double a = (*rows[row1])[col1];
    double b = (*rows[row2])[col1];
    double norm = sqrt(a*a + b*b);
    
    a = a/norm;
    b = b/norm;
    
    // @TODO make Givens a static? function    
    // update vector
    
    double en1 = (*c)[row1];
    double en2 = (*c)[row2];
    
    (*c)[row1] = a*en1 + b*en2;
    (*c)[row2] = -b*en1 + a*en2;  
    
    en1 = getEntry(row1, col1);
    en2 = getEntry(row2, col1);
    
    setEntry(row1, col1,  a*en1 + b*en2,true);
    dropFirst(row2);
    
    
    
    for(unsigned long j = col1 + 1; j <= rightIndex(row2); j++)
    {
        en1 = getEntry(row1, j);
        en2 = getEntry(row2, j);
        
        setEntry(row1, j,  a*en1 + b*en2,true);
        setEntry(row2, j, -b*en1 + a*en2,true);  
        
        //        cout<<"j "<<j<<" rowfill " <<rows[row1].getFill()<<endl;
    }
    
    
    for(int i = 0; i < rows[row1]->fillSize(); i++)
    {
        en1 = rows[row1]->getFill(i);
        en2 = rows[row2]->getFill(i);
        
        
        rows[row1]->setFill(i, a*en1 + b*en2);
        rows[row2]->setFill(i, -b*en1 + a*en2);  
        
    }
}


double FilledBandedOperator::rowDot(unsigned long row,unsigned long colm,unsigned long colM,vector<double> *r)
{
    double ret = 0;
    
    FilledRow *fr = rows[row];
    
    
    
    
    for(unsigned long j = colm; j <=colM; j++)
    {
        ret+=(*fr)[j]*(*r)[j];
    }
    
    return ret;
}
