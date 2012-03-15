//
//  main.cpp
//  Lecture 9
//
//  Created by Nis Sarup on 15/03/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "nr3.h"
//#include "roots.h"

template <class T> Doub rtnewt(T &funcd, const Doub x1, const Doub x2, const Doub xacc);
template <class T> Bool zbrac(T &func, Doub &x1, Doub &x2);

struct F{
    Doub operator()(Doub x) {
        return (x - cos(x));
    }
    
    double df(double x) {
        return (1 + sin(x));
    }
}f;

int main (int argc, const char * argv[])
{

    // insert code here...
    //cout << zbrac(f, -1000, 1000) << endl;
    double x1 = -1000.;
    double x2 = 1000.;
    zbrac(f, x1, x2);
    rtnewt<F>(f, x1, x2, 0.000001);
    std::cout << "Done!" << endl;
    return 0;
}

template <class T>
Doub rtnewt(T &funcd, const Doub x1, const Doub x2, const Doub xacc) {
    cout << setiosflags(ios::fixed) << setprecision(10);
    cout << "k|\t\tx(k)\t\t\td(k)\t\t\t\td/d" << endl;
    cout << "--------------------------------------------------------" << endl;
	const Int JMAX=20;
	Doub rtn=0.5*(x1+x2);
    Doub oldrtn;
    Doub oldd, d;
	for (Int j=0;j<JMAX;j++) {
		Doub f=funcd(rtn);
		Doub df=funcd.df(rtn);
		Doub dx=f/df;
		rtn -= dx;
        
        
        cout << j << "|\t\t" << rtn;
        if (j > 0) {
            d = rtn - oldrtn;
            cout << "\t" << d;
        }
        
        if (j > 1) {
            cout << "\t\t" << d / oldd;
        }
        
        cout << endl;
        oldrtn = rtn;
        oldd = d;
		if ((x1-rtn)*(rtn-x2) < 0.0)
			throw("Jumped out of brackets in rtnewt");
		if (abs(dx) < xacc) return rtn;
	}
	throw("Maximum number of iterations exceeded in rtnewt");
}

template <class T>
Bool zbrac(T &func, Doub &x1, Doub &x2)
{
	const Int NTRY=50;
	const Doub FACTOR=1.6;
	if (x1 == x2) throw("Bad initial range in zbrac");
	Doub f1=func(x1);
	Doub f2=func(x2);
	for (Int j=0;j<NTRY;j++) {
		if (f1*f2 < 0.0) return true;
		if (abs(f1) < abs(f2))
			f1=func(x1 += FACTOR*(x1-x2));
		else
			f2=func(x2 += FACTOR*(x2-x1));
	}
	return false;
}