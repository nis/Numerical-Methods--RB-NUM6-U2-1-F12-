import numpy
import scipy
import math

from scipy import linalg, matrix

def opgaveb1():
    A = numpy.matrix([[1.,4.,1.],[2.,1.,3.],[1.,1.,0.]]).transpose()
    x = numpy.matrix([9.,60.,7.]).transpose()
    
    alpha = numpy.linalg.solve(A,x)
    print alpha
    
def opgaveb2():
    A = numpy.matrix([  [math.sqrt(3./4.),              1./4.,                          -math.sqrt(3./2.)   ], 
                        [(3. - 2. * math.sqrt(3.))/8.,  (6. + math.sqrt(3.))/8.,        1./4.               ], 
                        [(2. + 3. * math.sqrt(3.))/8.,  (3. - 2. * math.sqrt(3.))/8.,   math.sqrt(3./4.)    ]
                    ]).transpose()
    A2 = numpy.matrix([  [math.sqrt(3.)/4.,              1./4.,                          -math.sqrt(3.)/2.  ], 
                        [(3. - 2. * math.sqrt(3.))/8.,  (6. + math.sqrt(3.))/8.,        1./4.               ], 
                        [(2. + 3. * math.sqrt(3.))/8.,  (3. - 2. * math.sqrt(3.))/8.,   math.sqrt(3.)/4.    ]
                    ]).transpose()
    
    x = numpy.matrix([8.,8.,8.]).transpose()
    
    R = A2.transpose()
    print R * x
    

#opgaveb1()
opgaveb2()