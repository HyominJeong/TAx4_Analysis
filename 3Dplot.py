import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# Create a sphere
r = 1
pi = np.pi
cos = np.cos
sin = np.sin
phi, theta = np.mgrid[0.0:pi:100j, 0.0:2.0*pi:100j]
x = r*sin(phi)*cos(theta)
y = r*sin(phi)*sin(theta)
z = r*cos(phi)

#Import data
data = np.genfromtxt('leb.txt')
r, theta, phi = np.hsplit(data, 3)
print r, theta, phi
theta = theta * pi / 180.0
phi = phi * pi / 180.0
xx = sin(phi)*cos(theta)
yy = sin(phi)*sin(theta)
zz = cos(phi)

#Set colours and render
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(
    x, y, z,  rstride=1, cstride=1, color='c', alpha=0.3, linewidth=0)

ax.scatter(xx,yy,zz,color="k",s=20)

ax.scatter(1.,0.,0.,color="r",s=20)
ax.scatter(0.,1.,0.,color="b",s=20)

ax.set_xlim([-1,1])
ax.set_ylim([-1,1])
ax.set_zlim([-1,1])
ax.set_aspect("equal", anchor="C")
plt.tight_layout()
plt.show()
