import funcHM
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import matplotlib.pyplot as plt
import sys

CR_energy = ([\
	5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5])
CR_ra = (np.array([\
	15,25,25,35,35,35,45,45,45,45,55,55,55,55,55,66,66,66,66,66,66]) +60+60\
	)* np.pi / 180
CR_dec = (np.array([\
	5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5])\
	)* np.pi / 180

#CR_events = [CR_energy, CR_ra, CR_dec]

#NofEvts, Ra_bin_edge, Dec_bin_edge = funcHM.nOfEvts_RaDec(CR_events)

#print NofEvts.transpose()

# CR flux = N of Evt / dOhm * Area * Time


# Draw
funcHM.drawFlux(CR_ra, CR_dec, binsize_deg=10)

