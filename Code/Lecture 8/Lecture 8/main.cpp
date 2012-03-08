//
//  main.cpp
//  Lecture 8
//
//  Created by Nis Sarup on 08/03/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include "nr3.h"
#include <iostream>
#include <ctime>    // For time()
#include <cstdlib>  // For srand() and rand()
#include <math.h>
#include <iomanip>

void do_pontius();
void printmat(MatDoub_I m);
MatDoub transpose(const MatDoub & x);
VecDoub matvecmult(MatDoub &A,VecDoub &b);
template<class T> NRmatrix<T> product(const NRmatrix<T>& A, const NRmatrix<T>& B);
template<class T> NRvector<T> product(const NRmatrix<T>& A, const NRvector<T>& b);
double length(const VecDoub b);
VecDoub divide_by_norm(const VecDoub b);
double dot(const VecDoub a, const VecDoub b);
double avg(const VecDoub x);
VecDoub scale(const VecDoub x);
VecDoub sum(const VecDoub a, const VecDoub b);
VecDoub difference(const VecDoub a, const VecDoub b);

int main (int argc, const char * argv[])
{

    do_pontius();
    std::cout << "Done!\n";
    return 0;
}

void do_pontius() {
    //float ys[40] = {.11019F, .21956F, .32949F, .43899F, .54803F, .65694F, .76562F, .87487F, .98292F, 1.09146F, 1.20001F, 1.30822F, 1.41599F, 1.52399F, 1.63194F, 1.73947F, 1.84646F, 1.95392F, 2.06128F, 2.16844F, .11052F, .22018F, .32939F, .43886F, .54798F, .65739F, .76596F, .87474F, .98300F, 1.09150F, 1.20004F, 1.30818F, 1.41613F, 1.52408F, 1.63159F, 1.73965F, 1.84696F, 1.95445F, 2.06177F, 2.16829F};
    double xs[40] = {150000, 300000, 450000, 600000, 750000, 900000, 1050000, 1200000, 1350000, 1500000, 1650000, 1800000, 1950000, 2100000, 2250000, 2400000, 2550000, 2700000, 2850000, 3000000, 150000, 300000, 450000, 600000, 750000, 900000, 1050000, 1200000, 1350000, 1500000, 1650000, 1800000, 1950000, 2100000, 2250000, 2400000, 2550000, 2700000, 2850000, 3000000};
    //
    MatDoub a(40,3);
    
    VecDoub x2(40);
    for (int i = 0; i < 40; i++) {
        x2[i] = xs[i];
    }
    
    VecDoub x3 = scale(x2);
    
    // Lav A-matricen
    for (int i = 0; i < 40; i++) {
        a[i][0] = 1;
        a[i][1] = x3[i];
        a[i][2] = (x3[i])*(x3[i]);
    }
    
    MatDoub ata = product(transpose(a), a);
    
    VecDoub x(3); // Initialize random number generator.
    x[0] = rand();
    x[1] = rand();
    x[2] = rand();
    
    double l = 0;
    VecDoub xold(3);
    
    VecDoub y = product(ata, x);
    
    cout << setiosflags(ios::fixed) << setprecision(8);
    cout << "k\t" << "l\t\t\t\t" << "d\t\t\t" << "d2\t\t\t" << endl;
    cout << "-----------------------------------------------" << endl;
    
    for (int i = 0; i < 10; i++) {
        xold = x;
        x = divide_by_norm(y);
        y = product(ata, x);
        l = dot(x, y);
        if (i == 0) {
            cout << i << "\t" << endl;
        } else if (i == 1) {
            cout << i << "\t";
            cout << l << "\t";
            cout << abs(length(difference(x, xold))) << "\t";
            cout << difference(x, xold)[2];
            cout << endl;
        } else if (i > 1) {
            cout << i << "\t";
            cout << l << "\t";
            cout << abs(length(difference(x, xold))) << "\t";
            cout << difference(x, xold)[2];
            cout << endl;
        }
    }
    
    //cout << l << endl;
}

VecDoub sum(const VecDoub a, const VecDoub b) {
    int n(a.size());
    if(n != b.size())
        throw "Vector sizes must be equal!";
    
    VecDoub result(n, 0.0);
    for(int i(0); i < n; i++)
        result[i] = a[i] + b[i];
    
    return result;
}

VecDoub difference(const VecDoub a, const VecDoub b) {
    int n(a.size());
    if(n != b.size())
        throw "Vector sizes must be equal!";
    
    VecDoub result(n, 0.0);
    for(int i(0); i < n; i++)
        result[i] = a[i] - b[i];
    
    return result;
}

VecDoub scale(const VecDoub x) {
    VecDoub res(x.size());
    double a = avg(x);
    
    for (int i = 0; i < x.size(); i++) {
        res[i] = x[i]/a;
    }
    return res;
}

double avg(const VecDoub x) {
    double res = 0;
    for (int i = 0; i < x.size(); i++) {
        res += x[i];
    }
    return res/x.size();
}

double dot(const VecDoub a, const VecDoub b) {
    double sum = 0;
    
    for (int i = 0; i < a.size(); i++) {
        sum += a[i] * b[i];
    }
    
    return sum;
}

VecDoub divide_by_norm(const VecDoub b) {
    double n = length(b);
    VecDoub res(b.size());
    for (int i = 0; i < b.size(); i++) {
        res[i] = 1/n * b[i];
    }
    return res;
}

double length(const VecDoub b) {
    double sum = 0;
    double result;
    
    for (int i = 0; i < b.size(); i++) {
        sum += b[i] * b[i];
    }
    
    result = sqrt(sum);
    return result;
}

// Matrix multiplication
template<class T> NRmatrix<T> product(const NRmatrix<T>& A, const NRmatrix<T>& B) {
    int r1(A.nrows()), r2(B.nrows()), c2(B.ncols());
    if(r2 != B.nrows())
        throw "Number of columns of first argument must equal number of rows of second argument!";
    
    NRmatrix<T> C(r1, c2, (T)0.0);
    for(int i = 0; i < r1; i++)
        for(int j = 0; j < c2; j++)
            for(int k = 0; k < r2; k++)
                C[i][j] += A[i][k] * B[k][j];
    
    return C;
}

// Multiplication between matrix and vector
template<class T> NRvector<T> product(const NRmatrix<T>& A, const NRvector<T>& b) {
    int r(A.nrows()), c(A.ncols());
    if(c != b.size())
        throw "Number of matrix columns must equal vector size!";
    
    NRvector<T> v(r, (T)0.0);
    for(int i = 0; i < r; i++)
        for(int j = 0; j < c; j++)
            v[i] += A[i][j] * b[j];
    
    return v;
}

VecDoub matvecmult(MatDoub &A,VecDoub &b)
{
    if(!A.ncols()==b.size()){cout<<"Wrong sizes in multiplication!"<<endl;return b;}
    VecDoub c(A.nrows());
    long double temp;
    for(int j=0;j<A.nrows();j++){
        temp=0.0;
        for(int k=0;k<A.ncols();k++){
            temp+=A[j][k]*b[k];
        }
        c[j]=temp;
    }
    return c;
}

MatDoub transpose(const MatDoub & x)
{
    MatDoub ret(x.ncols(), x.nrows());
    for (int i = 0; i < ret.nrows(); i++) {
        for (int j = 0; j < ret.ncols(); j++) {
            ret[i][j] = x[j][i];
        }
    }
    return ret;
}

void printmat(MatDoub_I m)
{
    for (int i = 0; i < m.nrows(); i++) {
        for (int j = 0; j < m.ncols(); j++) {
            cout << "\t" << setw(3) << m[i][j];
        }
        cout << endl;
    }
    cout << endl;
}