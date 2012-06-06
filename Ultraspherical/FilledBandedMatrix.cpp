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


FilledRow::FilledRow(int ind)
{
    index = ind;
}


int FilledRow::size()
{
    return entries.size();
}

double FilledRow::operator[](int j)
{
    if (FILLROWQ(j))
        return fill;
    else if(SHIFTROW(j) < 0)
        return 0;
    else
        return entries[SHIFTROW(j)];
}


void FilledRow::setEntry(int j,double val)
{
    setEntry(j, val,false);
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
    push_back(fill);
}


// increases Size so that row[k] is not a fill row
void FilledRow::increaseSize(int k)
{
    for(int i  = size(); i <= SHIFTROW(k); i++)
        increaseSize();
}


int FilledRow::leftIndex()
{
    return index;
}

int FilledRow::rightIndex()
{
    return index + size() - 1;
}




FilledBandedMatrix::FilledBandedMatrix(int low, int up)
{
    lowerIndex = low;
//    upper = up;
    increaseSize();
}





int FilledBandedMatrix::lower()
{
// Returns zero if the band starts on the diagonal
// -1 if there is one subdiagonal...
    return lowerIndex;
}
//int FilledBandedMatrix::upper()
//{
//    // Returns zero if the band ends on the diagonal
//    // 1 if there is one superdiagonal...
//    return rowSize - index - 1;
//}
//
//void FilledBandedMatrix::increaseUpper()
//{
//    vector<FilledRow>::iterator it;
//    for(it = rows.begin(); it < rows.end(); it++)
//    {
//        it->increaseSize();
//    }
//    
//    rowSize++;
//}

//void FilledBandedMatrix::increaseUpper(int k)
//{
//    for(int j = k; j >0; j--)
//    {
//        increaseUpper();
//    }
//}

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


FilledRow FilledBandedMatrix::createRow(int k)
{    
    double w = 15;
    
    if (k == 0) {
        FilledRow newrow(0);        
        newrow.setFill(1);
        return newrow;
    }
    else if (k == 1) {
        FilledRow newrow(0);                
        newrow.push_back(w);
        newrow.push_back(1);
        newrow.push_back(-.5*w);
        newrow.setFill(0);
        
        return newrow;
    } else {
        FilledRow newrow(k - 1);                
        newrow.push_back(.5*w);
        newrow.push_back(k);
        newrow.push_back(-.5*w);  
        newrow.setFill(0);        
        
        return newrow;        
    }    
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
    
    rows[i].setEntry(j,val);
}


FilledRow FilledBandedMatrix::back()
{
    return rows.back();
}


void FilledBandedMatrix::print()
{

    for(int i = 0; i < size(); i++)
    {
        for(int j = 0; j < columnSize(); j++)
        {
            cout << rows[i][j];
            cout << " ";
        }
        cout << endl;
    }
}





void FilledBandedMatrix::applyGivens(int row1, int row2, vector<double> *c)
{
    if(leftIndex(row2) < leftIndex(row1))
        throw "Givens should be applied left to right and top to bottom, in order";    
    
    
    int col1 = leftIndex(row1);
    
    double a = rows[row1][col1];
    double b = rows[row2][col1];
    double norm = sqrt(a*a + b*b);
    
    a = a/norm;
    b = b/norm;
    
// TODO make Givens a static? function    
// update vector
    
    double en1 = (*c)[row1];
    double en2 = (*c)[row2];
    
    (*c)[row1] = a*en1 + b*en2;
    (*c)[row2] = -b*en1 + a*en2;    
    

// update first row, deleting zeroed column    
    
    en1 = rows[row1][col1];
    en2 = rows[row2][col1];
    
    double val1 = a*en1+b*en2;
    double val2 = -b*en1+a*en2;            
    
    setEntry(row1, j, val1);
    setEntry(row2, j, val2);      
    
    
    
    for(int j = col1 + 1; j <= rightIndex(row2); j++)
    {
        en1 = rows[row1][j];
        en2 = rows[row2][j];
        
        double val1 = a*en1+b*en2;
        double val2 = -b*en1+a*en2;            
        
        setEntry(row1, j, val1);
        setEntry(row2, j, val2);  
        
//        cout<<"j "<<j<<" rowfill " <<rows[row1].getFill()<<endl;
    }
    
    double en1 = rows[row1][rightIndex(row1)+1];
    double en2 = rows[row2][rightIndex(row1)+1]; 

    
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
        ret+=fr[j]*(*r)[j];
    }
    
    return ret;
}