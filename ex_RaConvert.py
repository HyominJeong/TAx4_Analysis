from funcHM import RaConvt
import numpy as np

print "0 to pi -> -pi to pi"
a = np.linspace(0,2*np.pi,10)

print a*180/np.pi
b = RaConvt(a)
print b*180/np.pi

print "-pi to pi -> -pi to pi"
c = np.linspace(-np.pi, np.pi, 10)
print c*180/np.pi
d = RaConvt(c)
print d*180/np.pi