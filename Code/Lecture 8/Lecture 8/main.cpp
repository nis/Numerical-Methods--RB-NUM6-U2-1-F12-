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
#include "cholesky.h"
#include "ludcmp.h"
#include "svd.h"

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
void do_filip();
void potens_metode(MatDoub a, int precicion, int k);
void invers_potens_metode(MatDoub a, int precicion, int k);

int main (int argc, const char * argv[])
{

    do_pontius();
    do_filip();
    std::cout << "Done!\n";
    return 0;
}

void invers_potens_metode(MatDoub a, int precicion, int k) {
    VecDoub x(a.ncols());
    for (int i = 0; i < x.size(); i++) {
        x[i] = rand();
    }
    
    double l = 0;
    VecDoub xold(x.size());
    double d, dold, d2, d2old;
    
    //MatDoub ai(a.ncols(), a.ncols());
    //Cholesky c(a);
    //c.inverse(ai);
    //LUdcmp alu(a);
    //alu.inverse(ai);
    
    VecDoub y(x.size());
    
    //y = ai * x;
    SVD s(a);
    //s.solve(y, x);
    
    cout << setiosflags(ios::fixed) << setprecision(precicion);
    cout << "k\t" << "l\t\t\t\t" << "d\t\t\t\t" << "d2\t\t\t\t" << "d/d\t\t\t\t" << "d2/d2\t\t\t" << endl;
    cout << "--------------------------------------------------------------------------------" << endl;
    
    for (int i = 0; i < k; i++) {
        xold = x;
        s.solve(xold, y);
        x = divide_by_norm(y);
        l = dot(x, xold) / y.length();
        if (i == 0) {
            cout << i << "\t" << endl;
        } else if (i == 1) {
            d = abs(length(difference(x, xold)));
            d2 = difference(x, xold)[2];
            cout << i << "\t";
            cout << l << "\t";
            cout << d << "\t";
            cout << d2 << "\t";
            cout << endl;
        } else if (i > 1) {
            dold = d;
            d = abs(length(difference(x, xold)));
            d2old = d2;
            d2 = difference(x, xold)[2];
            cout << i << "\t";
            cout << l << "\t";
            cout << d << "\t";
            cout << d2 << "\t";
            cout << d/dold << "\t";
            cout << d2/d2old << "\t";
            cout << endl;
        }
    }
    cout << endl;
}

void potens_metode(MatDoub a, int precicion, int k) {
    VecDoub x(a.ncols());
    for (int i = 0; i < x.size(); i++) {
        x[i] = rand();
    }
    
    double l = 0;
    VecDoub xold(x.size());
    double d, dold, d2, d2old;
    
    VecDoub y = a * x;
    
    cout << setiosflags(ios::fixed) << setprecision(precicion);
    cout << "k\t" << "l\t\t\t\t" << "d\t\t\t\t" << "d2\t\t\t\t" << "d/d\t\t\t\t" << "d2/d2\t\t\t" << endl;
    cout << "--------------------------------------------------------------------------------" << endl;
    
    for (int i = 0; i < k; i++) {
        xold = x;
        x = divide_by_norm(y);
        y = a * x;
        l = dot(x, y);
        if (i == 0) {
            cout << i << "\t" << endl;
        } else if (i == 1) {
            d = abs(length(difference(x, xold)));
            d2 = difference(x, xold)[2];
            cout << i << "\t";
            cout << l << "\t";
            cout << d << "\t";
            cout << d2 << "\t";
            cout << endl;
        } else if (i > 1) {
            dold = d;
            d = abs(length(difference(x, xold)));
            d2old = d2;
            d2 = difference(x, xold)[2];
            cout << i << "\t";
            cout << l << "\t";
            cout << d << "\t";
            cout << d2 << "\t";
            cout << d/dold << "\t";
            cout << d2/d2old << "\t";
            cout << endl;
        }
    }
    cout << endl;
}

void do_filip() {
    double xs[82] = {-6.860120914, -4.324130045, -4.358625055, -4.358426747, -6.955852379, -6.661145254, -6.355462942, -6.118102026, -7.115148017, -6.815308569, -6.519993057, -6.204119983, -5.853871964, -6.109523091, -5.79832982, -5.482672118, -5.171791386, -4.851705903, -4.517126416, -4.143573228, -3.709075441, -3.499489089, -6.300769497, -5.953504836, -5.642065153, -5.031376979, -4.680685696, -4.329846955, -3.928486195, -8.56735134, -8.363211311, -8.107682739, -7.823908741, -7.522878745, -7.218819279, -6.920818754, -6.628932138, -6.323946875, -5.991399828, -8.781464495, -8.663140179, -8.473531488, -8.247337057, -7.971428747, -7.676129393, -7.352812702, -7.072065318, -6.774174009, -6.478861916, -6.159517513, -6.835647144, -6.53165267, -6.224098421, -5.910094889, -5.598599459, -5.290645224, -4.974284616, -4.64454848, -4.290560426, -3.885055584, -3.408378962, -3.13200249, -8.726767166, -8.66695597, -8.511026475, -8.165388579, -7.886056648, -7.588043762, -7.283412422, -6.995678626, -6.691862621, -6.392544977, -6.067374056, -6.684029655, -6.378719832, -6.065855188, -5.752272167, -5.132414673, -4.811352704, -4.098269308, -3.66174277, -3.2644011};
    
    VecDoub x(82);
    for (int i = 0; i < 82; i++) {
        x[i] = xs[i];
    }
    x = scale(x);
    
    MatDoub a(x.size(), 11);
    
    // Lav A-matricen
    for (int i = 0; i < x.size(); i++) {
        for (int j = 0; j < 11; j++) {
            a[i][j] = pow(x[i], (double)j);
        }
    }
    
    MatDoub ata = transpose(a) * a;
    
    cout << "Filip: Potens metode" << endl;
    potens_metode(ata, 18, 16);
    cout << "Filip: Invers potens metode" << endl;
    invers_potens_metode(ata, 18, 10);
}

void do_pontius() {
    double xs[40] = {150000, 300000, 450000, 600000, 750000, 900000, 1050000, 1200000, 1350000, 1500000, 1650000, 1800000, 1950000, 2100000, 2250000, 2400000, 2550000, 2700000, 2850000, 3000000, 150000, 300000, 450000, 600000, 750000, 900000, 1050000, 1200000, 1350000, 1500000, 1650000, 1800000, 1950000, 2100000, 2250000, 2400000, 2550000, 2700000, 2850000, 3000000};
    
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
    
    cout << "Pontius: Potens metode" << endl;
    potens_metode(ata, 10, 10);
    cout << "Pontius: Invers potens metode" << endl;
    invers_potens_metode(ata, 10, 10);
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