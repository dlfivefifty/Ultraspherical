//
//  BandedMatrix.cpp
//  Ultraspherical
//
//  Created by Sheehan Olver on 04/06/2012.
//  Copyright (c) 2012 School of Mathematics and Statistics, The University of Sydney. All rights reserved.
//

#include "FilledBandedMatrix.h"


#include <math.h>




FilledBandedMatrix::FilledBandedMatrix(int low,int up)
{
    index = -low;
    rowSize = up-low+1;
    
    increaseSize();
}

int FilledRow::size()
{
    return entries.size();
}

double FilledRow::operator[](int j)
{
    if (j >= size())
        return fill;
    else
        return entries[j];
}

void FilledRow::setEntry(int j,double val)
{    
    if (j >= size())
        setFill(val);
    else
        entries[j] = val;
}

void FilledRow::setFill(double val)
{
    fill = val;
}

double FilledRow::getFill()
{
    return fill;
}

void FilledRow::push_back(double val)
{    
    entries.push_back(val);
}

void FilledRow::increaseSize()
{
    entries.push_back(fill);
}




int FilledBandedMatrix::lower()
{
// Returns zero if the band starts on the diagonal
// -1 if there is one subdiagonal...
    return -index;
}
int FilledBandedMatrix::upper()
{
    // Returns zero if the band ends on the diagonal
    // 1 if there is one superdiagonal...
    return rowSize - index - 1;
}

void FilledBandedMatrix::increaseUpper()
{
    vector<FilledRow>::iterator it;
    for(it = rows.begin(); it < rows.end(); it++)
    {
        it->increaseSize();
    }
    
    rowSize++;
}

void FilledBandedMatrix::increaseUpper(int k)
{
    for(int j = k; j >0; j--)
    {
        increaseUpper();
    }
}

int FilledBandedMatrix::size()
{
    return rows.size();
}

int FilledBandedMatrix::columnSize()
{
    return size() + rows.back().size() + index;
}

int FilledBandedMatrix::columnSize(int col)
{
    return col - lower() + 1;
}

int FilledBandedMatrix::leftIndex(int row)
{
    return max(0,row+lower());
}

int FilledBandedMatrix::rightIndex(int row)
{
    return row+upper();
}


FilledRow FilledBandedMatrix::operator[] (int i)
{
    return rows[i];
}

double FilledBandedMatrix::getEntry(int i,int j)
{    
    if(i >= size()) 
        return 0;
    else if (j - i + index < 0)
        return 0;
    else
        return rows[i][j - i + index];
    
}

FilledRow FilledBandedMatrix::createRow(int k)
{
    FilledRow newrow;
    
    double w = 15;
    
    if (k == 0) {
        for(int i = 0; i < rowSize; i++) {
            newrow.push_back(1);
        }
        newrow.setFill(1);
    }
    else if (k == 1) {
        newrow.push_back(w);
        newrow.push_back(1);
        newrow.push_back(-.5*w);
        for(int i = 3; i < rowSize; i++) {  
            newrow.push_back(0);
        }
        newrow.setFill(0);
    } else {
        newrow.push_back(.5*w);
        newrow.push_back(k);
        newrow.push_back(-.5*w);  
        for(int i = 3; i < rowSize; i++) {  
            newrow.push_back(0);
        }        
        newrow.setFill(0);        
    }

    
    
    return newrow;
    
}

void FilledBandedMatrix::increaseSize()
{
    rows.push_back(createRow(size()));
}


void FilledBandedMatrix::setEntry(int i,int j,double val)
{    
    for(int k = size(); k <= i; k++)
    {
        rows.push_back(createRow(k));
    }
    
    rows[i].setEntry(j - i + index,val);    
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


void FilledBandedMatrix::applyGivens(int row1, int row2, double a, double b)
{
//    cout <<"rightind: "<<rightIndex(row1)<<endl;
    
    for(int j = leftIndex(row2); j <= rightIndex(row1); j++)
    {
        double en1 = getEntry(row1, j);
        double en2 = getEntry(row2, j);
        
        double val1 = a*en1+b*en2;
        double val2 = -b*en1+a*en2;            
        
        setEntry(row1, j, val1);
        setEntry(row2, j, val2);  
        
//        cout<<"j "<<j<<" rowfill " <<rows[row1].getFill()<<endl;
    }
    
    double en1 = getEntry(row1, rightIndex(row1)+1);
    double en2 = getEntry(row2, rightIndex(row1)+1); 

    
    double val1 = a*en1+b*en2;
    double val2 = -b*en1+a*en2;            

    
    for(int j = rightIndex(row1)+1; j <= rightIndex(row2)+1; j++)
    {        
        setEntry(row1, j, val1);
        setEntry(row2, j, val2);         
    }
}


double FilledBandedMatrix::rowDot(int row,int colm,int colM,vector<double> *r)
{
    double ret = 0;
    
    FilledRow fr = rows[row];
    
    
    
    
    for(int j = colm; j <=colM; j++)
    {
        ret+=fr[j - row + index]*(*r)[j];
    }
    
    return ret;
}