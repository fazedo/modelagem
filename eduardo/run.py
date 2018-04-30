import time
import calor
import calor2
import calor3
import sys
import numpy
 
#a = numpy.matrix(
#    [[1, -2, 0, 0, 0],
#    [-2, 1, -2, 0, 0],
#    [0, -2, 1, -2, 0],
#    [0, 0, -2, 1, -2],
#    [0, 0, 0, -2, 1]]
#)
#b = numpy.array([1.0, 1, 1, 1, 1])
#d1 = numpy.array([-2.0, -2, -2, -2])
#d3 = numpy.array([-2.0, -2, -2, -2])
#d2 = numpy.array([1.0, 1, 1, 1, 1])
#c = numpy.linalg.solve(a, b)
#print(c)
 
#c2 = calor3.TDMAsolver(d1, d2, d3, b)
#print(c2)
#print(d1, d2, d3, b)
#sys.exit(0)
 
start_time = time.time()
[x, u] = calor2.calor(1.0, 0.1, 201, 2001)
elapsed_time = time.time() - start_time
print('elapsed_time:', elapsed_time)
#print(u)
 
start_time = time.time()
calor.calor(1.0, 0.1, 201, 2001)
elapsed_time = time.time() - start_time
 
print('elapsed_time:', elapsed_time)
#print(x)
#print(u)
 
start_time = time.time()
[x, u] = calor3.calor(1.0, 0.1, 201, 2001)
elapsed_time = time.time() - start_time
print('elapsed_time:', elapsed_time)
print(u)
