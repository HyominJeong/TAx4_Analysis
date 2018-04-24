import numpy as np

delta_deg = 10.

def dohm(theta_deg, phi_deg):
	result = delta_deg * np.pi / 180 * (-np.cos((theta_deg+0.5*delta_deg)*np.pi/180) + np.cos((theta_deg-0.5*delta_deg)*np.pi/180))
	return result

#print dohm(85.,5.)
dohms = [dohm(5+10*n, 5) for n in range(9)]
print sum(dohms) * 36
