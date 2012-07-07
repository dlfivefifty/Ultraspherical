
#include "MultiplicationMatrices.h"
#include <math.h>
#include "AdaptiveQR.h"
// TODO: Better to construct the Hankel vector in reverse so that it can be pop_backed as the ConvertMult vector is formed.


FilledRow zeroRowGenerator(int k)
{
    FilledRow zero(0,RowFiller::dirichlet(0, 0));
    zero.increaseSize();
    return zero;
}



double applyConversion(vector<double> row,int shift,int k)
{
    int ind = shift;
    
//    printvec(row);
    
    double ret = 0;

    if(ind >= (int)row.size())
        return ret;
    
    if( ind >= 0)
    {
        if(k == 0)
            ret = row[ind];
        else
            ret = row[ind]*.5/(k+1);
    }
    
    ind+=2;
    
    if(ind >= (int)row.size())
        return ret;
    
    if ( ind >= 0)
        ret += -row[ind]*(.5/(k+3.) + .5/(k+1));
    
    ind+=2;
    
    if(ind >= (int)row.size())
        return ret;
    
    if(ind >= 0)
        ret += row[ind]*.5/(k+3);
    
    return ret; 
}

double applyConversion(vector<double> cin,int k)
{
    return applyConversion(cin,0, k);
}

double zeroFirstApplyConversion(vector<double> cin,int k)
{
//    printvec(cin);
    cin.insert(cin.begin(),0);
//    printvec(cin);    
    return applyConversion(cin,0, k);
}





ToeplitzMatrix::ToeplitzMatrix(vector<double> a) : FilledBandedMatrix(a.size() - 1, NULL)
{
    vector<double> halved;
    
    vector<double>::iterator it = a.begin();
    
    halved.push_back(*it);
    for (++it; it < a.end(); it++)
    {
        halved.push_back(.5*(*it));
    }
    
    
    
    push_back(FilledRow(0,halved,RowFiller::dirichlet(0,0)));
    
    it = halved.begin();
//    it2++;
    
    vector<double> halvedadd = halved;
    
    

    for(++it; it < halved.end(); it++)
    {   
        halvedadd.insert(halvedadd.begin(),*it);
        push_back(FilledRow(0,halvedadd,RowFiller::dirichlet(0,0)));
    }   
    
    
    rowEntries = halvedadd;
}

FilledRow ToeplitzMatrix::createRow(int k)
{
    return FilledRow(k-(rowEntries.size()-1)/2,rowEntries,RowFiller::dirichlet(0, 0));
}


MultiplicationMatrix::MultiplicationMatrix(vector<double> a) : FilledBandedMatrix(a.size() - 1, NULL)
{
    vector<double> halved;
    
    vector<double>::iterator it = a.begin();
    
    halved.push_back(*it);
    for (++it; it < a.end(); it++)
    {
        halved.push_back(.5*(*it));
    }
    
    
    
    push_back(FilledRow(0,halved,RowFiller::dirichlet(0,0)));
    
    it = halved.begin();
    //    it2++;
    
    vector<double> halvedadd = halved;
    
    
    
    for(++it; it < halved.end(); it++)
    {   
        halvedadd.insert(halvedadd.begin(),*it);
        
        vector<double> newrow;
        
        vector<double>::iterator hankelit = it;
        
        for (vector<double>::iterator init = halvedadd.begin(); init < halvedadd.end(); ++init) {
            if (hankelit < halved.end()) 
            {
                newrow.push_back((*init) + (*hankelit));
                hankelit++;
            }
            else
                newrow.push_back((*init));
        }
        
        
        push_back(FilledRow(0,newrow,RowFiller::dirichlet(0,0)));
    }   
    
    
    rowEntries = halvedadd;
}

FilledRow MultiplicationMatrix::createRow(int k)
{
    return FilledRow(k-(rowEntries.size()-1)/2,rowEntries,RowFiller::dirichlet(0, 0));
}



ConvertToeplitzMatrix::ConvertToeplitzMatrix(vector<double> a) : FilledBandedMatrix(a.size() - 1, NULL)
{
    vector<double> halved;
    
    vector<double>::iterator it = a.begin();
    
    halved.push_back(*it);
    for (++it; it < a.end(); it++)
    {
        halved.push_back(.5*(*it));
    }
    
    
    
//    push_back(FilledRow(0,halved,RowFiller::dirichlet(0,0)));
    
    it = halved.begin();
    //    it2++;
    
    vector<double> halvedadd = halved;
    
    
    vector<double> firstrow;
    
    
    firstrow.push_back(applyConversion(halved, 0));
    
    
    for(++it; it < halved.end(); it++)
    {   
        halvedadd.insert(halvedadd.begin(),*it);
        firstrow.push_back(applyConversion(halvedadd, 0));        
    }
    
    vector<double> halvedaddwithzeros = halvedadd;
    for(int i = 0; i < 4; i++)
    {
        halvedaddwithzeros.insert(halvedaddwithzeros.begin(), 0);
        firstrow.push_back(applyConversion(halvedaddwithzeros, 0));                
    }
    
    
    push_back(FilledRow(0,firstrow,RowFiller::dirichlet(0,0)));    
    
    

    
    
//    for (int strt = 1; strt < a.size(); strt++) {
//        vector<double> newrow;
//        
//        for(int i = a.size() + strt + 3; i >=0; i--)
//        {
//            newrow.push_back(applyConversion(halvedaddwithzeros,i, strt));
//        }
//        
//        push_back(FilledRow(0,newrow,RowFiller::dirichlet(0,0)));    
//    }
//    
    for(int i = 1; i < halved.size(); i++)
    {
        vector<double> newrow;        
        for (int j =0; j < halved.size() + i + 4; j++) {
            newrow.push_back(applyConversion(halvedaddwithzeros, halved.size() +i + 3 -j,i));            
        }
        push_back(FilledRow(0,newrow,RowFiller::dirichlet(0,0)));    
    }    
    
    
    rowEntries = halvedadd;
}

FilledRow ConvertToeplitzMatrix::createRow(int k)
{
    

    vector<double> newrow;
    
    for (int i = rowEntries.size()-1; i >= -4; i--) {
        newrow.push_back(applyConversion(rowEntries,i, k)); 
    }
    
    return FilledRow(k-(rowEntries.size()-1)/2,newrow,RowFiller::dirichlet(0, 0));
}



ConvertHankelMatrix::ConvertHankelMatrix(vector<double> a) : FilledBandedMatrix(a.size() - 1, NULL)
{
    vector<double> halved;
    
    vector<double>::iterator it = a.begin();
    
    halved.push_back(*it);
    for (++it; it < a.end(); it++)
    {
        halved.push_back(.5*(*it));
    }
    
    
    
    //    push_back(FilledRow(0,halved,RowFiller::dirichlet(0,0)));
    
//    it = halved.begin();
    //    it2++;
    
    vector<double> halvedremove = halved;
    
    vector<double> firstrow;
    
    while(halvedremove.size() >0)
    {
        halvedremove.erase(halvedremove.begin());    
        firstrow.push_back(zeroFirstApplyConversion(halvedremove, 0));
    }
    
    push_back(FilledRow(0,firstrow,RowFiller::dirichlet(0,0)));    
    
    
    for(int i = 2; i < halved.size(); i++)
    {
        vector<double> newrow;        
        for (int j =i-1; j < halved.size(); j++) {
            newrow.push_back(applyConversion(halved, j,i-1));            
        }
        push_back(FilledRow(0,newrow,RowFiller::dirichlet(0,0)));    
    }
    

}


FilledRow ConvertHankelMatrix::createRow(int k)
{
    
    
    vector<double> newrow;
    newrow.push_back(0);
    
    return FilledRow(1,newrow,RowFiller::dirichlet(0, 0));
}



ConvertMultiplicationMatrix::ConvertMultiplicationMatrix(vector<double> a) : FilledBandedMatrix(a.size() - 1, NULL)
{
    vector<double> halved;
    
    vector<double>::iterator it = a.begin();
    
    halved.push_back(*it);
    for (++it; it < a.end(); it++)
    {
        halved.push_back(.5*(*it));
    }
    
    
    
    //    push_back(FilledRow(0,halved,RowFiller::dirichlet(0,0)));
    
    it = halved.begin();
    //    it2++;
    
    vector<double> halvedadd = halved;
    
    
    vector<double> firstrow;
    
    
    vector<double> halvedremove = halved;
    
//    while(halvedremove.size() >0)
//    {
//        halvedremove.erase(halvedremove.begin());    
//        firstrow.push_back(zeroFirstApplyConversion(halvedremove, 0));
//    }

    
    halvedremove.erase(halvedremove.begin());            
    firstrow.push_back(applyConversion(halved, 0)+zeroFirstApplyConversion(halvedremove, 0));
    
    
    for(++it; it < halved.end(); it++)
    {   
        halvedremove.erase(halvedremove.begin()); 
        
        halvedadd.insert(halvedadd.begin(),*it);
        firstrow.push_back(applyConversion(halvedadd, 0)+zeroFirstApplyConversion(halvedremove, 0));        
    }
    
    vector<double> halvedaddwithzeros = halvedadd;
    for(int i = 0; i < 4; i++)
    {
        halvedaddwithzeros.insert(halvedaddwithzeros.begin(), 0);
        firstrow.push_back(applyConversion(halvedaddwithzeros, 0));                
    }
    
    
    push_back(FilledRow(0,firstrow,RowFiller::dirichlet(0,0)));    
    
    
    
    
    
//    for (int strt = 1; strt < a.size(); strt++) {
//        vector<double> newrow;
//        
//        for(int i = a.size() + strt + 3; i >=0; i--)
//        {
//            newrow.push_back(applyConversion(halvedaddwithzeros,i, strt));
//        }
//        
//        push_back(FilledRow(0,newrow,RowFiller::dirichlet(0,0)));    
//    }//i + 3 + a.size(), i + 2 + a.size(),   a.size() - i
    
    for(int i = 1; i < halved.size(); i++)
    {
        vector<double> newrow;        
        for (int j =0; j < halved.size() - i; j++) {
            newrow.push_back(applyConversion(halved, j + i,i) + applyConversion(halvedaddwithzeros, halved.size() +i + 3 -j,i));            
        }
        
        for (int j =halved.size() - i; j < halved.size() + i + 4; j++) {
            newrow.push_back(applyConversion(halvedaddwithzeros, halved.size() +i + 3 -j,i));            
        }        
        
        
        push_back(FilledRow(0,newrow,RowFiller::dirichlet(0,0)));    
    }
        
    
    rowEntries = halvedadd;
}



FilledRow ConvertMultiplicationMatrix::createRow(int k)
{
    
    
    vector<double> newrow;
    
    for (int i = rowEntries.size()-1; i >= -4; i--) {
        newrow.push_back(applyConversion(rowEntries,i, k)); 
    }
    
    return FilledRow(k-(rowEntries.size()-1)/2,newrow,RowFiller::dirichlet(0, 0));
}





DirichletD2ConvertMultiplicationMatrix::DirichletD2ConvertMultiplicationMatrix(vector<double> a) : FilledBandedMatrix(a.size() - 1, NULL)
{
    vector<double> halved;
    
    vector<double>::iterator it = a.begin();
    
    
    push_back(FilledRow(0,RowFiller::dirichlet(1,0)));    
    push_back(FilledRow(0,RowFiller::dirichlet(0,1)));        
    
    halved.push_back(*it);
    for (++it; it < a.end(); it++)
    {
        halved.push_back(.5*(*it));
    }
    
    

    
    vector<double> halvedadd = halved;
    
    
    vector<double> firstrow;
    
    
    vector<double> halvedremove = halved;
    
    //    while(halvedremove.size() >0)
    //    {
    //        halvedremove.erase(halvedremove.begin());    
    //        firstrow.push_back(zeroFirstApplyConversion(halvedremove, 0));
    //    }
    
    
    halvedremove.erase(halvedremove.begin());            
    firstrow.push_back(applyConversion(halved, 0)+zeroFirstApplyConversion(halvedremove, 0));
    
    
    for(int i = 1; i < halved.size(); i++)
    {   
        halvedremove.erase(halvedremove.begin()); 
        
        halvedadd.insert(halvedadd.begin(),halved[i]);
        if(i == 2)
            firstrow.push_back(applyConversion(halvedadd, 0)+zeroFirstApplyConversion(halvedremove, 0) + 4);        
        else
            firstrow.push_back(applyConversion(halvedadd, 0)+zeroFirstApplyConversion(halvedremove, 0));                    
    }
    
    vector<double> halvedaddwithzeros = halvedadd;
    for(int i = 0; i < 4; i++)
    {
        halvedaddwithzeros.insert(halvedaddwithzeros.begin(), 0);
        firstrow.push_back(applyConversion(halvedaddwithzeros, 0));                
    }
    
    
    push_back(FilledRow(0,firstrow,RowFiller::dirichlet(0,0)));    
    
    
    
    
    
    //    for (int strt = 1; strt < a.size(); strt++) {
    //        vector<double> newrow;
    //        
    //        for(int i = a.size() + strt + 3; i >=0; i--)
    //        {
    //            newrow.push_back(applyConversion(halvedaddwithzeros,i, strt));
    //        }
    //        
    //        push_back(FilledRow(0,newrow,RowFiller::dirichlet(0,0)));    
    //    }//i + 3 + a.size(), i + 2 + a.size(),   a.size() - i
    
    for(int i = 1; i < halved.size(); i++)
    {
        vector<double> newrow;        
        for (int j =0; j < halved.size() - i; j++) {
            if(j == 2 + i)
                newrow.push_back(applyConversion(halved, j + i,i) + applyConversion(halvedaddwithzeros, halved.size() +i + 3 -j,i) + 4 + 2*i);        
            else
                newrow.push_back(applyConversion(halved, j + i,i) + applyConversion(halvedaddwithzeros, halved.size() +i + 3 -j,i));                              
        }
        
        for (int j =halved.size() - i; j < halved.size() + i + 4; j++) {            
            if(j == 2 + i)
                newrow.push_back(applyConversion(halvedaddwithzeros, halved.size() +i + 3 -j,i) + 4 + 2*i);        
            else
                newrow.push_back(applyConversion(halvedaddwithzeros, halved.size() +i + 3 -j,i));                
        }        
        
        
        push_back(FilledRow(0,newrow,RowFiller::dirichlet(0,0)));    
    }
    
    
    rowEntries = halvedadd;
}



FilledRow DirichletD2ConvertMultiplicationMatrix::createRow(int k)
{
    
    
    vector<double> newrow;
    
//    for (int i = rowEntries.size()+3; i >= 0; i--) {
//        newrow.push_back(applyConversion(rowEntries,i-4, k-2)); 
//    }
    
    for (int i = 0; i <= rowEntries.size()+3; i++) {
        if(i + k-2-(rowEntries.size()-1)/2 == k)
            newrow.push_back(applyConversion(rowEntries,rowEntries.size()+3-i-4, k-2) + 4 + 2*(k-2)); 
        else
            newrow.push_back(applyConversion(rowEntries,rowEntries.size()+3-i-4, k-2));             
    }    
    
    return FilledRow(k-2-(rowEntries.size()-1)/2,newrow,RowFiller::dirichlet(0, 0));
}







FilledBandedMatrix hankelOperator(vector<double> a)
{
    FilledBandedMatrix ret(a.size(),&zeroRowGenerator);
    
    ret.increaseSize();
    
    for (vector<double>::iterator it = a.begin(); it != a.end(); ++it) {
        FilledRow row(0,RowFiller::dirichlet(0, 0));        
        for (vector<double>::iterator init = it; init != a.end(); ++init) {
            row.push_back(.5*(*init));
        };
        
        
        ret.push_back(row);
    };
    
    
    return ret;
}




// S_0 * M_0[a]. 
vector<double> ConvertMult(vector<double> a, int k,int n)
{
    a[0] = 2*a[0]; 
    vector<double> hankel; //hankel diagonal
    vector<double> mdiag; //diagonal
    hankel = ConvertHankel(a,k);
    
    double dg;
    int s,t;  
    s = k*(k >=0) + -k*(k<0); // starting index, abs(k). 
    t = (s-2)*(s-2>=0) + (2-s)*(s-2<0); // abs(s-2); 
    for (int j = 0 ; j < s + 3 ; ++j)
        a.push_back(0); // append some zeros. 
    
    dg = ((a[s] - a[s+2])/2)*(k<=0) + ((a[s] - a[t])/2)*(k>0) ;
    
    //first element is special case.
    mdiag.push_back((dg + (k>=0)*a[s]/2 + hankel[0])/2); 
    
    // calculate Toeplitz part of conversion * multiplication. 
    for (vector<double>::iterator it=hankel.begin()+1; it<hankel.end(); ++it)
    {
        mdiag.push_back((dg+*it)/2);
    };
    
    //hankel part over, constant on diagonal. 
    for (int j = mdiag.size(); j < n ; ++j)
    {
        mdiag.push_back(dg/2);
    };
    
    return mdiag; 
}

// constructs transformed hankel matrix.

vector<double> ConvertHankel(vector<double> a,int k) //returns the kth diagonal.
{
    
    vector<double> mdiag;
    int s; 
    s = k*(k >=0) + -k*(k<0); // starting index, abs(k). 
    a.push_back(0); a.push_back(0); 
    //first element is special case.
    mdiag.push_back((a[s]*(k<0) - a[s+2])/2); 
    
    // calculate hankel part of conversion * multiplication. 
    for (vector<double>::iterator it = a.begin()+s+2; it < a.end()-2; it = it+2) 
    {
        mdiag.push_back((*it - *(it+2))/2);
    };
    
    return mdiag; 
}

// M_0[a]*D_1
//vector<double> MultDiff(vector<double> a, int k,int n)
//{

//return
//}

/*
 // M_1[a]*D_1
 vector<double> MultDiff(vector<double> a, int k,int n)
 {
 a[0]=2*a[0]; 
 vector<double> mdiag; 
 vector<double> hankel; //hankel diagonal
 int m; //starting multiple
 int s; //abs(k)
 
 
 hankel = HankelDiff(a,k); 
 m = (k<=0) + (k+1)*(k>0); 
 s = k*(k>=0) + -k*(k<0); 
 
 mdiag.push_back(0);
 for (vector<double>::iterator it = a.begin()+s+1; it < a.end(); it = it+2) 
 {
 mdiag.push_back(m*(*it));
 m=m+1;
 };
 
 return mdiag;
 }
 */

// Returns the kth diagonal of Hankel * diff matrix. 
vector<double> HankelDiff(vector<double> a, int k)
{
    vector<double> sum;
    vector<double> mdiag; 
    double tot=0; 
    int s; //abs(k)
    int len = a.size(); 
    
    //m = (k<=0) + (k+1)*(k>0); 
    s = k*(k>=0) + -k*(k<0); 
    int sgn = (s+1) % 2; 
    
    // make a of odd or even length depending on k. 
    int j=0; 
    while ( ((len+j) % 2 + s) % 2 ) 
    {
        a.push_back(0); 
        j=j+1;
    }
    
    // even/odd sum. 
    for (vector<double>::iterator it = a.begin()+sgn; it<a.end(); it=it+2) 
    {
        tot = tot + (*it);
        sum.push_back(tot);  // with diff multipliers.
    };
    
    // apply diff multipliers and flip.
    
    double mult = k*(k>0); 
    for (vector<double>::iterator it = sum.begin()+s-(k>1); it<sum.end()-sgn*(k>0); ++it) 
    {
        mdiag.push_back(mult*(*it));  // with diff multipliers.
        mult = mult + 1; 
    };
    return mdiag;
}

