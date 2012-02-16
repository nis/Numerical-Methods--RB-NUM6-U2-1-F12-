//
//  main.cpp
//  Lecture 3
//
//  Created by Nis Sarup on 16/02/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "nr3.h"
#include "ludcmp.h"
#include "svd.h"
#include "math.h"

void do_lu(MatDoub a, VecDoub b);
void do_svd(MatDoub a, VecDoub b, double t);
void printmat(MatDoub_I m);
MatDoub transpose(const MatDoub & x);
VecDoub matvecmult(MatDoub &A,VecDoub &b);
template<class T> NRmatrix<T> product(const NRmatrix<T>& A, const NRmatrix<T>& B);
template<class T> NRvector<T> product(const NRmatrix<T>& A, const NRvector<T>& b);
double sse(MatDoub a, VecDoub b, VecDoub x);

int main (int argc, const char * argv[])
{
    int m = 82;
    double ys[82] = {0.8116F, 0.9072F, 0.9052F, 0.9039F, 0.8053F, 0.8377F, 0.8667F, 0.8809F, 0.7975F, 0.8162F, 0.8515F, 0.8766F, 0.8885F, 0.8859F, 0.8959F, 0.8913F, 0.8959F, 0.8971F, 0.9021F, 0.909F, 0.9139F, 0.9199F, 0.8692F, 0.8872F, 0.89F, 0.891F, 0.8977F, 0.9035F, 0.9078F, 0.7675F, 0.7705F, 0.7713F, 0.7736F, 0.7775F, 0.7841F, 0.7971F, 0.8329F, 0.8641F, 0.8804F, 0.7668F, 0.7633F, 0.7678F, 0.7697F, 0.77F, 0.7749F, 0.7796F, 0.7897F, 0.8131F, 0.8498F, 0.8741F, 0.8061F, 0.846F, 0.8751F, 0.8856F, 0.8919F, 0.8934F, 0.894F, 0.8957F, 0.9047F, 0.9129F, 0.9209F, 0.9219F, 0.7739F, 0.7681F, 0.7665F, 0.7703F, 0.7702F, 0.7761F, 0.7809F, 0.7961F, 0.8253F, 0.8602F, 0.8809F, 0.8301F, 0.8664F, 0.8834F, 0.8898F, 0.8964F, 0.8963F, 0.9074F, 0.9119F, 0.9228};
    double xs[82] = {-6.860120914, -4.324130045, -4.358625055, -4.358426747, -6.955852379, -6.661145254, -6.355462942, -6.118102026, -7.115148017, -6.815308569, -6.519993057, -6.204119983, -5.853871964, -6.109523091, -5.79832982, -5.482672118, -5.171791386, -4.851705903, -4.517126416, -4.143573228, -3.709075441, -3.499489089, -6.300769497, -5.953504836, -5.642065153, -5.031376979, -4.680685696, -4.329846955, -3.928486195, -8.56735134, -8.363211311, -8.107682739, -7.823908741, -7.522878745, -7.218819279, -6.920818754, -6.628932138, -6.323946875, -5.991399828, -8.781464495, -8.663140179, -8.473531488, -8.247337057, -7.971428747, -7.676129393, -7.352812702, -7.072065318, -6.774174009, -6.478861916, -6.159517513, -6.835647144, -6.53165267, -6.224098421, -5.910094889, -5.598599459, -5.290645224, -4.974284616, -4.64454848, -4.290560426, -3.885055584, -3.408378962, -3.13200249, -8.726767166, -8.66695597, -8.511026475, -8.165388579, -7.886056648, -7.588043762, -7.283412422, -6.995678626, -6.691862621, -6.392544977, -6.067374056, -6.684029655, -6.378719832, -6.065855188, -5.752272167, -5.132414673, -4.811352704, -4.098269308, -3.66174277, -3.2644011};
    
    MatDoub a(m,11);
    VecDoub b(m);
    
    // Lav A-matricen
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < 11; j++) {
            a[i][j] = pow(xs[i], (double)j);
        }
    }
    
    // Lav b-matricen
    for (int i = 0; i < m; i++) {
        b[i] = ys[i];
    }
    
    do_lu(a, b);
    
    do_svd(a, b, -1);
    
    do_svd(a, b, 0.0000000001);
    
    std::cout << "Done!" <<endl;
    return 0;
}

void do_svd(MatDoub a, VecDoub b, double t) {
    
    cout << "Solving with SVD:" << endl;
    
    // Lav resultat-vektor:
    VecDoub x(11);
    
    SVD svd(a);
    svd.solve(b, x, t);
    
    // Print out the awesome results!
    for (int i = 0; i < 11; i++) {
        cout << "x" << i << ": \t" << x[i] << endl;
    }
    
    cout << endl << "SSE:\t" << sse(a, b, x) << endl;
    
    cout << endl;
}

void do_lu(MatDoub a, VecDoub b) {
    
    cout << "Solving with LU:" << endl;
    
    MatDoub at(11, 82);
    
    // Transponér A
    at = transpose(a);
    
    // Lav A-stjerne: AT·A
    MatDoub astjerne = product(at, a);
    
    // Lav b-stjerne: AT·b
    VecDoub bstjerne = product(at, b);
    
    
    // Lav resultat-vektor:
    VecDoub x(11);
    
    LUdcmp alu(astjerne);
    alu.solve(bstjerne, x);
    
    // Print out the awesome results!
    for (int i = 0; i < 11; i++) {
        cout << "x" << i << ": \t" << x[i] << endl;
    }
    
    cout << endl << "SSE:\t" << sse(a, b, x) << endl;
    cout << endl;
}

double sse(MatDoub a, VecDoub b, VecDoub x) {
    double sse = 0;
    VecDoub ssetemp = product(a, x);
    for (int i = 0; i < ssetemp.size(); i++) {
        sse += pow(ssetemp[i] - b[i], 2); 
    }
    return sse;
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