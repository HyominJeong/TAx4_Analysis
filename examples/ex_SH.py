from scipy.special import sph_harm
import numpy as np
import matplotlib.pyplot as plt

#theta = np.linspace(0, 2*np.pi, 100)
theta = np.pi
phi = np.linspace(0, np.pi, 100)
m = -1
n = 1
Ymn = sph_harm(m,n,theta, phi)
Ym11 = 0.5 * np.sqrt(3./2/np.pi) * np.exp(-1.j * theta) * np.sin(phi)
plt.plot(phi, Ymn, 'r-')
plt.plot(phi, Ym11, 'b:')
plt.show()