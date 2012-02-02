//
//  main.cpp
//  Lecture_1
//
//  Created by Nis Sarup on 02/02/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "nr3.h"
#include "ludcmp.h"

int main (int argc, const char * argv[])
{
    const Int n = 3;
    MatDoub a(n,n);
    VecDoub b(n),x(n);
    
    /*
     A-matrix:
     2  4   8
     4  8   4
     -2 10  5
     
     b-matrix:
     1
     2
     3
    */
    
    // Entering data
    a[0][0] = 2;
    a[0][1] = 4;
    a[0][2] = 8;
    a[1][0] = 4;
    a[1][1] = 8;
    a[1][2] = 4;
    a[2][0] = -2;
    a[2][1] = 10;
    a[2][2] = 5;
    
    b[0] = 1;
    b[1] = 2;
    b[2] = 3;
    
    
    LUdcmp alu(a);
    
    alu.solve(b, x);
    
    
    // Outputting results
    std::cout << "x0 = " << x[0] << endl;
    std::cout << "x1 = " << x[1] << endl;
    std::cout << "x2 = " << x[2] << endl;
    
    
    std::cout << "[ " << alu.lu[0][0] << " " << alu.lu[0][1] << " "  << alu.lu[0][2] << " "  << "]" << endl;
    std::cout << "[ " << alu.lu[1][0] << " " << alu.lu[1][1] << " "  << alu.lu[1][2] << " "  << "]" << endl;
    std::cout << "[ " << alu.lu[2][0] << " " << alu.lu[2][1] << " "  << alu.lu[2][2] << " "  << "]" << endl;

    std::cout << "Done!\n";
    return 0;
}

