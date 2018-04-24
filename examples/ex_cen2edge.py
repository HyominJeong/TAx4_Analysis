import funcHM
import numpy as np

a = np.zeros((3,3), dtype=np.dtype(object))
a[0,0] = np.array([0.0, 0.0])
a[0,1] = np.array([0.0, 1.0])
a[0,2] = np.array([0.0, 2.0])
a[1,0] = np.array([1.0, 0.0])
a[1,1] = np.array([1.0, 1.0])
a[1,2] = np.array([1.0, 2.0])
a[2,0] = np.array([2.0, 0.0])
a[2,1] = np.array([2.0, 1.0])
a[2,2] = np.array([2.0, 2.0])

print "Input array is \n", a

result = funcHM.bin_edge2cent(a)

print result
print len(result[0])