import numpy
import scipy

from scipy import linalg, matrix

vectora = [[1,-3], [2,-5]]
vectorb = [[1,2], [2,5], [3,8]]
vectorc = [[1,1,1], [1,0,0], [0,1,0]]
vectord = [[1,2,4], [1,1,1], [2,3,6]]
vectore = [[1,2,4], [1,3,5], [-1,-1,-3]]
vectorf = [[2,4,6,8], [1,1,1,1]]
vectorg = [[1,3,6,9],[2,4,6,8],[5,13,24,35]]
vectorh = [[1,1,1,1],[2,2,2,3],[3,3,4,5],[1,1,0,0],[0,0,1,1]]

def rank(A, eps=1e-12):
    u, s, vh = numpy.linalg.svd(A)
    return len([x for x in s if abs(x) > eps])
    
def is_liniearly_independent(a):
    m = numpy.matrix(a)
    m = m.transpose()
    if (len(a) > rank(m)):
        return False
    else:
        return True
        
def null(A, eps=1e-15):
    u, s, vh = scipy.linalg.svd(A)
    null_mask = (s <= eps)
    null_space = scipy.compress(null_mask, vh, axis=0)
    return scipy.transpose(null_space)

def opgave_1():
    print 'Opgave 1:'
    print '----------------------------'
    print('Vector a: ' + str(is_liniearly_independent(vectora)))
    print('Vector b: ' + str(is_liniearly_independent(vectorb)))
    print('Vector c: ' + str(is_liniearly_independent(vectorc)))
    print('Vector d: ' + str(is_liniearly_independent(vectord)))
    print('Vector e: ' + str(is_liniearly_independent(vectore)))
    print('Vector f: ' + str(is_liniearly_independent(vectorf)))
    print('Vector g: ' + str(is_liniearly_independent(vectorg)))
    print('Vector h: ' + str(is_liniearly_independent(vectorh)))
    print '----------------------------'
    
def opgave_2():
    print 'Opgave 2:'
    print '----------------------------'
    print 'A\'s rang: ' + str(rank(numpy.matrix(vectora).transpose()))
    print 'C\'s rang: ' + str(rank(numpy.matrix(vectorc).transpose()))
    print 'D\'s rang: ' + str(rank(numpy.matrix(vectord).transpose()))
    print 'E\'s rang: ' + str(rank(numpy.matrix(vectore).transpose()))
    print 'F\'s rang: ' + str(rank(numpy.matrix(vectorf).transpose()))
    print 'G\'s rang: ' + str(rank(numpy.matrix(vectorg).transpose()))
    print '----------------------------'

def opgave_3():
    print 'Opgave 3:'
    print '----------------------------'
    print 'B\'s rang: ' + str(rank(numpy.matrix(vectorb)))
    print 'H\'s rang: ' + str(rank(numpy.matrix(vectorh)))
    print '----------------------------'

def opgave_4():
    print 'Opgave 4:'
    print '----------------------------'
    print 'ATA\'s fra opgave 2'
    print 'A\'s rang: ' + str(rank(numpy.matrix(vectora).transpose() * numpy.matrix(vectora)))
    print 'C\'s rang: ' + str(rank(numpy.matrix(vectorc).transpose() * numpy.matrix(vectorc)))
    print 'D\'s rang: ' + str(rank(numpy.matrix(vectord).transpose() * numpy.matrix(vectord)))
    print 'E\'s rang: ' + str(rank(numpy.matrix(vectore).transpose() * numpy.matrix(vectore)))
    print 'F\'s rang: ' + str(rank(numpy.matrix(vectorf).transpose() * numpy.matrix(vectorf)))
    print 'G\'s rang: ' + str(rank(numpy.matrix(vectorg).transpose() * numpy.matrix(vectorg)))
    print ''
    print 'ATA\'s fra opgave 3'
    print 'B\'s rang: ' + str(rank(numpy.matrix(vectorb).transpose() * numpy.matrix(vectorb)))
    print 'H\'s rang: ' + str(rank(numpy.matrix(vectorh).transpose() * numpy.matrix(vectorh)))
    print '----------------------------'

def opgave_5():
    print 'Opgave 5:'
    print '----------------------------'
    #A = numpy.matrix([[1,0,0,0,2],[0,0,3,0,0],[0,0,0,0,0],[0,4,0,0,0]])
    A = numpy.matrix(vectorg)
    #print rank(A)
    print A
    print ''
    U, S, V = numpy.linalg.svd(A,full_matrices=0)
    print U
    print ''
    print S
    print ''
    print V
    print '----------------------------'


#opgave_1()
#print numpy.matrix(vectorg).transpose()
#opgave_2()
#opgave_3()
#print numpy.matrix(vectorh)
#opgave_4()
#
opgave_5()