import funcHM
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import matplotlib.pyplot as plt
import sys

############################
# Make Homogeneous CR evts #
############################

nOfSample = 5000

CRH_ra, CRH_dec = funcHM.CRHradec(nOfSample, 1)

##########################################
# Differential exposure from declination #
##########################################
# Not applied

CRD_ra, CRD_dec = CRH_ra, CRH_dec

###################################
# Get histogram of position delta #
###################################
nbin = 1000
hist, binEdge = funcHM.histPosDelta(CRD_ra, CRD_dec, nbin)

# Draw histogram
bincenters = 0.5*(binEdge[1:] + binEdge[:-1])
binsizes = 0.5*(binEdge[1:] - binEdge[:-1])
menStd = np.sqrt(hist)

plt.errorbar(bincenters, hist, menStd, binsizes)
plt.show()
