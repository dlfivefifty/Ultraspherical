
#include "MultiplicationMatrices.h"
#include <math.h>
#include "AdaptiveQR.h"
// TODO: Better to construct the Hankel vector in reverse so that it can be pop_backed as the ConvertMult vector is formed.


FilledRow zeroRowGenerator(unsigned long k)
{
    FilledRow zero(0,RowFiller::dirichlet(0, 0));
    zero.increaseSize();
    return zero;
}



double applyConversion(vector<double> *row,int shift,unsigned long k)
{
    int ind = shift;
    
//    printvec(row);
    
    double ret = 0;

    if(ind >= (int)row->size())
        return ret;
    
    if( ind >= 0)
    {
        if(k == 0)
            ret = (*row)[ind];
        else
            ret = (*row)[ind]*.5/(k+1);
    }
    
    ind+=2;
    
    if(ind >= (int)row->size())
        return ret;
    
    if ( ind >= 0)
        ret += -(*row)[ind]*(.5/(k+3.) + .5/(k+1));
    
    ind+=2;
    
    if(ind >= (int)row->size())
        return ret;
    
    if(ind >= 0)
        ret += (*row)[ind]*.5/(k+3);
    
    return ret; 
}

double applyConversion(vector<double> *cin,int k)
{
    return applyConversion(cin,0, k);
}

double zeroFirstApplyConversion(vector<double> *row)
{
    
    
////    printvec(cin);
//    cin.insert(cin.begin(),0);
////    printvec(cin);    
//    return applyConversion(cin,0, k);
//    
//  TODO Fix following    
//    int ind = 0; k = 0;
    
    //    printvec(row);
    
    double ret = 0;
    
//    if(ind >= (int)row->size())
//        return ret;
//    
//        if(k == 0)
//            ret = (*row)[ind];
//        else
//            ret = (*row)[ind]*.5/(k+1);
    
    int ind=1;
    
    if(ind >= (int)row->size())
        return ret;
    
    ret += -(*row)[ind]*(.5/(3.) + .5/(1));
    
    ind+=2;
    
    if(ind >= (int)row->size())
        return ret;
    
    
    ret += (*row)[ind]*.5/(3);
    
    return ret;     
}





//ToeplitzMatrix::ToeplitzMatrix(vector<double> a) : FilledBandedMatrix(a.size() - 1, NULL)
//{
//    vector<double> *halved = new vector<double>;
//    
//    vector<double>::iterator it = a.begin();
//    
//    halved->push_back(*it);
//    for (++it; it < a.end(); it++)
//    {
//        halved->push_back(.5*(*it));
//    }
//    
//    
//    
//    push_back(FilledRow(0,halved,RowFiller::dirichlet(0,0)));
//    
//    it = halved->begin();
////    it2++;
//    
//    vector<double> *halvedadd = new vector<double>(halved->begin(),halved->end());
//    
//
//    for(++it; it < halved->end(); it++)
//    {   
//        halvedadd->insert(halvedadd->begin(),*it);
//        push_back(FilledRow(0,halvedadd,RowFiller::dirichlet(0,0)));
//    }   
//    
//    
//    rowEntries = halvedadd;
//}
//
//FilledRow ToeplitzMatrix::createRow(int k)
//{
//    return FilledRow(k-(rowEntries->size()-1)/2,rowEntries,RowFiller::dirichlet(0, 0));
//}
//
//
//MultiplicationMatrix::MultiplicationMatrix(vector<double> a) : FilledBandedMatrix(a.size() - 1, NULL)
//{
//    vector<double> *halved = new vector<double>;
//    
//    vector<double>::iterator it = a.begin();
//    
//    halved->push_back(*it);
//    for (++it; it < a.end(); it++)
//    {
//        halved->push_back(.5*(*it));
//    }
//    
//    
//    
//    push_back(FilledRow(0,halved,RowFiller::dirichlet(0,0)));
//    
//    it = halved->begin();
//    //    it2++;
//    
//    vector<double> *halvedadd = new vector<double>(halved->begin(),halved->end());    
//    
//    
//    for(++it; it < halved->end(); it++)
//    {   
//        halvedadd->insert(halvedadd->begin(),*it);
//        
//        vector<double> newrow;
//        
//        vector<double>::iterator hankelit = it;
//        
//        for (vector<double>::iterator init = halvedadd.begin(); init < halvedadd->end(); ++init) {
//            if (hankelit < halved->end()) 
//            {
//                newrow.push_back((*init) + (*hankelit));
//                hankelit++;
//            }
//            else
//                newrow.push_back((*init));
//        }
//        
//        
//        push_back(FilledRow(0,&newrow,RowFiller::dirichlet(0,0)));
//    }   
//    
//    
//    rowEntries = halvedadd;
//}
//
//FilledRow MultiplicationMatrix::createRow(int k)
//{
//    return FilledRow(k-(rowEntries.size()-1)/2,rowEntries,RowFiller::dirichlet(0, 0));
//}
//
//
//
//ConvertToeplitzMatrix::ConvertToeplitzMatrix(vector<double> a) : FilledBandedMatrix(a.size() - 1, NULL)
//{
//    vector<double> halved;
//    
//    vector<double>::iterator it = a.begin();
//    
//    halved.push_back(*it);
//    for (++it; it < a.end(); it++)
//    {
//        halved.push_back(.5*(*it));
//    }
//    
//    
//    
////    push_back(FilledRow(0,halved,RowFiller::dirichlet(0,0)));
//    
//    it = halved.begin();
//    //    it2++;
//    
//    vector<double> halvedadd = halved;
//    
//    
//    vector<double> firstrow;
//    
//    
//    firstrow.push_back(applyConversion(halved, 0));
//    
//    
//    for(++it; it < halved.end(); it++)
//    {   
//        halvedadd.insert(halvedadd.begin(),*it);
//        firstrow.push_back(applyConversion(halvedadd, 0));        
//    }
//    
//    vector<double> halvedaddwithzeros = halvedadd;
//    for(int i = 0; i < 4; i++)
//    {
//        halvedaddwithzeros.insert(halvedaddwithzeros.begin(), 0);
//        firstrow.push_back(applyConversion(halvedaddwithzeros, 0));                
//    }
//    
//    
//    push_back(FilledRow(0,firstrow,RowFiller::dirichlet(0,0)));    
//    
//    
//
//    
//    
////    for (int strt = 1; strt < a.size(); strt++) {
////        vector<double> newrow;
////        
////        for(int i = a.size() + strt + 3; i >=0; i--)
////        {
////            newrow.push_back(applyConversion(halvedaddwithzeros,i, strt));
////        }
////        
////        push_back(FilledRow(0,newrow,RowFiller::dirichlet(0,0)));    
////    }
////    
//    for(int i = 1; i < halved.size(); i++)
//    {
//        vector<double> newrow;        
//        for (int j =0; j < halved.size() + i + 4; j++) {
//            newrow.push_back(applyConversion(halvedaddwithzeros, halved.size() +i + 3 -j,i));            
//        }
//        push_back(FilledRow(0,newrow,RowFiller::dirichlet(0,0)));    
//    }    
//    
//    
//    rowEntries = halvedadd;
//}
//
//FilledRow ConvertToeplitzMatrix::createRow(int k)
//{
//    
//
//    vector<double> newrow;
//    
//    for (int i = rowEntries.size()-1; i >= -4; i--) {
//        newrow.push_back(applyConversion(rowEntries,i, k)); 
//    }
//    
//    return FilledRow(k-(rowEntries.size()-1)/2,newrow,RowFiller::dirichlet(0, 0));
//}
//
//
//
//ConvertHankelMatrix::ConvertHankelMatrix(vector<double> a) : FilledBandedMatrix(a.size() - 1, NULL)
//{
//    vector<double> halved;
//    
//    vector<double>::iterator it = a.begin();
//    
//    halved.push_back(*it);
//    for (++it; it < a.end(); it++)
//    {
//        halved.push_back(.5*(*it));
//    }
//    
//    
//    
//    //    push_back(FilledRow(0,halved,RowFiller::dirichlet(0,0)));
//    
////    it = halved.begin();
//    //    it2++;
//    
//    vector<double> halvedremove = halved;
//    
//    vector<double> firstrow;
//    
//    while(halvedremove.size() >0)
//    {
//        halvedremove.erase(halvedremove.begin());    
//        firstrow.push_back(zeroFirstApplyConversion(halvedremove, 0));
//    }
//    
//    push_back(FilledRow(0,firstrow,RowFiller::dirichlet(0,0)));    
//    
//    
//    for(int i = 2; i < halved.size(); i++)
//    {
//        vector<double> newrow;        
//        for (int j =i-1; j < halved.size(); j++) {
//            newrow.push_back(applyConversion(halved, j,i-1));            
//        }
//        push_back(FilledRow(0,newrow,RowFiller::dirichlet(0,0)));    
//    }
//    
//
//}
//
//
//FilledRow ConvertHankelMatrix::createRow(int k)
//{
//    
//    
//    vector<double> newrow;
//    newrow.push_back(0);
//    
//    return FilledRow(1,newrow,RowFiller::dirichlet(0, 0));
//}
//
//
//
//ConvertMultiplicationMatrix::ConvertMultiplicationMatrix(vector<double> a) : FilledBandedMatrix(a.size() - 1, NULL)
//{
//    vector<double> halved;
//    
//    vector<double>::iterator it = a.begin();
//    
//    halved.push_back(*it);
//    for (++it; it < a.end(); it++)
//    {
//        halved.push_back(.5*(*it));
//    }
//    
//    
//    
//    //    push_back(FilledRow(0,halved,RowFiller::dirichlet(0,0)));
//    
//    it = halved.begin();
//    //    it2++;
//    
//    vector<double> halvedadd = halved;
//    
//    
//    vector<double> firstrow;
//    
//    
//    vector<double> halvedremove = halved;
//    
////    while(halvedremove.size() >0)
////    {
////        halvedremove.erase(halvedremove.begin());    
////        firstrow.push_back(zeroFirstApplyConversion(halvedremove, 0));
////    }
//
//    
//    halvedremove.erase(halvedremove.begin());            
//    firstrow.push_back(applyConversion(halved, 0)+zeroFirstApplyConversion(halvedremove, 0));
//    
//    
//    for(++it; it < halved.end(); it++)
//    {   
//        halvedremove.erase(halvedremove.begin()); 
//        
//        halvedadd.insert(halvedadd.begin(),*it);
//        firstrow.push_back(applyConversion(halvedadd, 0)+zeroFirstApplyConversion(halvedremove, 0));        
//    }
//    
//    vector<double> halvedaddwithzeros = halvedadd;
//    for(int i = 0; i < 4; i++)
//    {
//        halvedaddwithzeros.insert(halvedaddwithzeros.begin(), 0);
//        firstrow.push_back(applyConversion(halvedaddwithzeros, 0));                
//    }
//    
//    
//    push_back(FilledRow(0,firstrow,RowFiller::dirichlet(0,0)));    
//    
//    
//    
//    
//    
////    for (int strt = 1; strt < a.size(); strt++) {
////        vector<double> newrow;
////        
////        for(int i = a.size() + strt + 3; i >=0; i--)
////        {
////            newrow.push_back(applyConversion(halvedaddwithzeros,i, strt));
////        }
////        
////        push_back(FilledRow(0,newrow,RowFiller::dirichlet(0,0)));    
////    }//i + 3 + a.size(), i + 2 + a.size(),   a.size() - i
//    
//    for(int i = 1; i < halved.size(); i++)
//    {
//        vector<double> newrow;        
//        for (int j =0; j < halved.size() - i; j++) {
//            newrow.push_back(applyConversion(halved, j + i,i) + applyConversion(halvedaddwithzeros, halved.size() +i + 3 -j,i));            
//        }
//        
//        for (int j =halved.size() - i; j < halved.size() + i + 4; j++) {
//            newrow.push_back(applyConversion(halvedaddwithzeros, halved.size() +i + 3 -j,i));            
//        }        
//        
//        
//        push_back(FilledRow(0,newrow,RowFiller::dirichlet(0,0)));    
//    }
//        
//    
//    rowEntries = halvedadd;
//}
//
//
//
//FilledRow ConvertMultiplicationMatrix::createRow(int k)
//{
//    
//    
//    vector<double> newrow;
//    
//    for (int i = rowEntries.size()-1; i >= -4; i--) {
//        newrow.push_back(applyConversion(rowEntries,i, k)); 
//    }
//    
//    return FilledRow(k-(rowEntries.size()-1)/2,newrow,RowFiller::dirichlet(0, 0));
//}





FilledBandedMatrix *DirichletD2ConvertMultiplicationMatrix(vector<double> a) 
{
    PlusRowAdder *pl = new PlusRowAdder(new ConversionRowAdder());
    pl->push_back(new DerivativeRowAdder());
    FilledBandedMatrix *ret = new FilledBandedMatrix(-2,new ShiftRowAdder(pl,-2));
    
    FilledRow *drrow = new FilledRow(0,RowFiller::dirichlet(1,0));
    
    drrow->increaseSize();
    ret->push_back(drrow);
    
    drrow = new FilledRow(0,RowFiller::dirichlet(0,1));
    drrow->increaseSize();
    ret->push_back(drrow);
    
    return ret;
}








