//
//  main.cpp
//  Varmeplader
//
//  Created by Nis Sarup on 12/04/12.
//  Copyright (c) 2012 Sarup. All rights reserved.
//

#include <iostream>
#include "nr3.h"
#include "math.h"
#include "ludcmp.h"

template <class T> double int_trapez(T &f, int N, double a, double b);
template <class T> void sol(T &F, int N, double a, double b, double e1, double e2, double T1, double T2);

struct F{
    double operator()(double h, double x, double y, double e) {
        double res = 0;
        res = 1 - e;
        double top = pow(1.0, 2.0);
        double bottom = pow(pow(1.0, 2.0) + pow((x - y), 2.0), (3.0/2.0));
        res = res * 0.5 * (top / bottom);
        res = h * res;
        return res;
    }
}Pl;


int main (int argc, const char * argv[])
{

    // insert code here...
    int N = 4;
    double a = -0.5;
    double b = 0.5;
    double e1 = 0.80;
    double e2 = 0.60;
    double t1 = 1000.0;
    double t2 = 500.0;
    sol<F>(Pl, N, a, b, e1, e2, t1, t2);
    std::cout << "Done!" << endl;
    return 0;
}

template <class T> void sol(T &F, int N, double a, double b, double e1, double e2, double T1, double T2) {
    
    // Make the A matrix
    MatDoub A(N*2, N*2);
    for (int i = 0; i <= (N*2) - 1; i++) {
        A[i][i] = 1;
    }
    
    // Calculate h
    double h = (b - a)/N;
    
    // Make and populate the x and y vectors
    VecDoub x(N);
    VecDoub y(N);
    for (int i = 0; i < x.size(); i++) {
        x[i] = a + i * h;
        cout << x[i] << endl;
        y[i] = a + i * h;
    }
    
    
    // Do the first set
    for (int ypos = 0; ypos < N; ypos++) {
        for (int xpos = N; xpos < (N * 2); xpos++) {
            double ht = 0;
            if (xpos == N || xpos == (N * 2) - 1) {
                ht = h / 2;
            } else {
                ht = h;
            }
            A[ypos][xpos] = F(ht, x[xpos], y[ypos], e1);
        }
    }
    
    // Do the second set
    for (int ypos = N; ypos < (N * 2); ypos++) {
        for (int xpos = 0; xpos < N; xpos++) {
            double ht = 0;
            if (ypos == N || ypos == (N * 2) - 1) {
                ht = h / 2;
            } else {
                ht = h;
            }
            A[ypos][xpos] = F(ht, x[xpos], y[ypos], e2);
        }
    }
    //A.print();
    
    // Make the b vector
    VecDoub bvec(N * 2);
    for (int i = 0; i < bvec.size(); i++) {
        if (i < N) {
            bvec[i] = e1 * (1.712 * pow(10.0, -9.0)) * pow(T1, 4.0);
        } else {
            bvec[i] = e2 * (1.712 * pow(10.0, -9.0)) * pow(T2, 4.0);
        }
    }
    //bvec.print();
    
    VecDoub res(N * 2);
    
    LUdcmp alu(A);
    alu.solve(bvec, res);
    for (int i = 0; i < res.size(); i++) {
        cout << res[i] << endl;
    }
}

/**********************************************************************
 * Integral calculation using trapez method
 */
template <class T>
double int_trapez(T &f, int N, double a, double b){
    double sum = 0;
    double x = 0;
    double h = (b - a)/N;
    for (int i = 1; i < N; i++){
        x = a + i*h;
        sum += 2*f(x);
    }
    sum += f(a) + f(b);
    return h*sum/2;
}