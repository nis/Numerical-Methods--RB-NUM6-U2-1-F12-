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
    const Int n = 4;
    MatDoub a(n,n);
    VecDoub b(n),x(n);
    
    LUdcmp alu(a);
    
    // insert code here...
    std::cout << "Hello, World!\n";
    return 0;
}

