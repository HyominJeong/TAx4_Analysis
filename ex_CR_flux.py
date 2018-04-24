import funcHM
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import matplotlib.pyplot as plt
import sys

CRasciiFile = sys.argv[1]
MCasciiFile = sys.argv[2]
eThr = float(sys.argv[3])
eLim = float(sys.argv[4])
degbin = int(sys.argv[5])

# CR events
CR_yymmdd, CR_hhmmss, CR_energy, CR_ra, CR_dec = funcHM.read_TAascii(CRasciiFile, eThr, eLim, form='wiki')
#CR_yymmdd, CR_hhmmss, CR_energy, CR_ra, CR_dec = funcHM.read_TAascii(CRasciiFile, eThrPower, form='anisotropy')
'''
txt_radec = open("radec_from_wiki.txt",'w')
for i in range(30):
	print>>txt_radec, ("%6d\t%6.6f\t%.2f\t%.2f\t%.2f" % (CR_yymmdd[i], CR_hhmmss[i], CR_energy[i], CR_ra[i], CR_dec[i]))
'''
CR_ra = np.pi - CR_ra
CR_events = [CR_energy, CR_ra, CR_dec]



#NofEvts, Ra_bin_edge, Dec_bin_edge = funcHM.nOfEvts_RaDec(CR_events)

#print NofEvts.transpose()

# CR flux = N of Evt / dOhm * Area * Time


# Draw
flux, Pos_Ra, Pos_Dec = funcHM.drawFlux(CR_ra, CR_dec, binsize_deg=degbin)

print len(flux[0]), len(Pos_Ra[0]), len(Pos_Dec[0])
print flux[0]
print Pos_Ra[0]
print Pos_Dec[0]

# Print a_lm
l_arr = []
c_l_arr = []
for l in range(100):
	#l = 1
	c_l = 0
	for m in range(-l, l+1):
		a_lm = funcHM.a_lm(flux, Pos_Ra, Pos_Dec, l, m)
		c_l += np.absolute(a_lm) ** 2
		#print np.sqrt(a_lm * np.conj(a_lm))
	#print l, c_l
	l_arr.append(l)
	c_l_arr.append(c_l / (2*l +1))

plt.plot(l_arr, c_l_arr)
#plt.set_xscale("log")
plt.yscale("log")
plt.show()

plt.plot(l_arr, c_l_arr)
#plt.set_xscale("log")
#plt.yscale("log")
plt.show()