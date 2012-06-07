//
//  BandedMatrix.cpp
//  Ultraspherical
//
//  Created by Sheehan Olver on 04/06/2012.
//  Copyright (c) 2012 School of Mathematics and Statistics, The University of Sydney. All rights reserved.
//

#include "FilledBandedMatrix.h"


#include <math.h>

#define SHIFTROW(j)(j - index)
#define FILLROWQ(j)(j - index >= size())


double oneFiller(int k)
{
    return 1;
}

double alternatingFiller(int k)
{
    return 1-2*(k % 2);
}

RowFiller::RowFiller()
{
    fill = 0;
    fillGenerator = &oneFiller;
}

RowFiller::RowFiller(double fil, double (*filr)(int))   
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

double RowFiller::getEntry(int k)
{
   return  fill*fillGenerator(k);
}


double RowFiller::fillGenerate(int col)
{
    return fillGenerator(col);
}



FilledRow FilledRow::leftDirichlet()
{
    RowFiller flr(1,&alternatingFiller);
    vector<RowFiller> fillers;
    fillers.push_back(flr);
    FilledRow newrow(0,fillers);        
    return newrow;
}

FilledRow FilledRow::rightDirichlet()
{
    RowFiller flr(0,&oneFiller);
    vector<RowFiller> fillers;
    fillers.push_back(flr);
    FilledRow newrow(0,fillers);     
    return newrow;
}

FilledRow::FilledRow(int ind)
{
    index = ind;
}


FilledRow::FilledRow(int ind,vector<RowFiller> flr)
{
    index = ind;
    fillers = flr;    
}

FilledRow::FilledRow(int ind,vector<RowFiller> flr,int size)
{
    index = ind;
    fillers = flr;
    
    entries.resize(size);
}


int FilledRow::size()
{
    return entries.size();
}

double FilledRow::operator[](int j)
{
    //    cout<<"size(): "<<entries.size()<<endl;
    
    if (FILLROWQ(j)) {
        double ret = 0;
        for(vector<RowFiller>::iterator it = fillers.begin(); it < fillers.end(); it++)
            ret+=it->getEntry(j);
        
        return ret;
    }
    else if(SHIFTROW(j) < 0)
        return 0;
    else
        return entries[SHIFTROW(j)];
}


void FilledRow::setEntry(int j,double val)
{
    setEntry(j,val,false);
}
void FilledRow::setEntry(int j,double val,bool increasesize)
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
    
    entries[SHIFTROW(j)] = val;
}

void FilledRow::setFill(double val)
{
    fillers[0].setFill(val);
}

double FilledRow::getFill()
{
    return fillers[0].getFill();
}

double FilledRow::fillGenerate(int col)
{
    return fillers[0].fillGenerate(col);
}

void FilledRow::push_back(double val)
{    
    entries.push_back(val);
}

void FilledRow::increaseSize()
{
    push_back((*this)[size()]);
}

// increases Size so that row[k] is not a fill row
void FilledRow::increaseSize(int k)
{
    for(int i  = size(); i <= SHIFTROW(k); i++)
        increaseSize();
}

void FilledRow::dropFirst()
{
    if(size() == 0)
        throw "No first to drop";
    
    entries.erase(entries.begin());
    
    index++;
}


int FilledRow::leftIndex()
{
    return index;
}

int FilledRow::rightIndex()
{
    return index + size() - 1;
}



FilledBandedMatrix::FilledBandedMatrix(int low,FilledRow (*rgen)(int))
{
    lowerIndex = low;
    rowGenerator = rgen;
    
    increaseSize();
}




int FilledBandedMatrix::lower()
{
    // Returns zero if the band starts on the diagonal
    // -1 if there is one subdiagonal...
    return lowerIndex;
}


int FilledBandedMatrix::size()
{
    return rows.size();
}

int FilledBandedMatrix::columnSize()
{
    return rows.back().rightIndex();
}

int FilledBandedMatrix::columnSize(int col)
{
    return col - lower() + 1;
}

int FilledBandedMatrix::leftIndex(int row)
{
    return max(0,rows[row].leftIndex());
}

int FilledBandedMatrix::rightIndex(int row)
{
    return rows[row].rightIndex();
}


FilledRow FilledBandedMatrix::operator[] (int i)
{
    return rows[i];
}

double FilledBandedMatrix::getEntry(int i,int j)
{    
    if(i >= size()) 
        return 0;
    else
        return rows[i][j];
    
}

FilledRow FilledBandedMatrix::createRow(int k)
{
    return (*rowGenerator)(k);
}

void FilledBandedMatrix::increaseSize()
{
    rows.push_back(createRow(size()));
}

void FilledBandedMatrix::dropFirst(int row)
{
    rows[row].dropFirst();
}


void FilledBandedMatrix::setEntry(int i,int j,double val)
{
    setEntry(i,j,false);
}

void FilledBandedMatrix::setEntry(int i,int j,double val,bool increasesize)
{    
    for(int k = size(); k <= i; k++)
    {
        rows.push_back(createRow(k));
    }
    
    rows[i].setEntry(j,val,increasesize);    
}



void FilledBandedMatrix::print()
{
    
    for(int i = 0; i < size(); i++)
    {
        for(int j = 0; j < columnSize(); j++)
        {
            cout << getEntry(i,j);
            cout << " ";
        }
        cout << endl;
    }
}


void FilledBandedMatrix::applyGivens(int row1, int row2, vector<double> *c)
{
    int col1 = row1;    
    if(leftIndex(row2) < leftIndex(row1))
        throw "Givens should be applied left to right and top to bottom, in order";    
    if (leftIndex(row2) > col1)
        throw "Cannot change elements to the left";
    if (rightIndex(row1) > rightIndex(row2))
        throw "Probably a bug: rightIndex(row1) > rightIndex(row2)";
    if (row2 > size())
        throw "Probably a bug: row2 > size";
    
    
    double a = (rows[row1])[col1];
    double b = (rows[row2])[col1];
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
    
    
    
    for(int j = col1 + 1; j <= rightIndex(row2); j++)
    {
        en1 = getEntry(row1, j);
        en2 = getEntry(row2, j);
        
        setEntry(row1, j,  a*en1 + b*en2,true);
        setEntry(row2, j, -b*en1 + a*en2,true);  
        
        //        cout<<"j "<<j<<" rowfill " <<rows[row1].getFill()<<endl;
    }
    
    en1 = rows[row1].getFill();
    en2 = rows[row2].getFill();
    
    rows[row1].setFill( a*en1 + b*en2);
    rows[row2].setFill(-b*en1 + a*en2);  
}


double FilledBandedMatrix::rowDot(int row,int colm,int colM,vector<double> *r)
{
    double ret = 0;
    
    FilledRow fr = rows[row];
    
    
    
    
    for(int j = colm; j <=colM; j++)
    {
        ret+=fr[j]*(*r)[j];
    }
    
    return ret;
}