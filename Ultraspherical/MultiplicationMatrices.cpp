
#include "MultiplicationMatrices.h"
#include <math.h>

// TODO: Better to construct the Hankel vector in reverse so that it can be pop_backed as the ConvertMult vector is formed.



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

