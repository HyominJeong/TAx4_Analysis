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
flux, Pos_Ra, Pos_Dec = funcHM.drawFlux(CR_ra, CR_dec, binsize_deg=15)

print len(flux[0]), len(Pos_Ra[0]), len(Pos_Dec[0])
print flux[0]
print Pos_Ra[0]
print Pos_Dec[0]

