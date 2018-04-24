import funcHM
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import matplotlib.pyplot as plt
import sys

##########################################################
# Check arrival direction gap distribution               #
# by Hyomin Jeong                                        #
# Last modify: 2018.04.17.                               #
##########################################################

# How to use:
# python CR_gap_dist.py CR_event(ascii file) MC_event(ascii file) Minimum_energy(float, EeV unit)


CRasciiFile = sys.argv[1]
MCasciiFile = sys.argv[2]
eThrPower = float(sys.argv[3])

# CR events
CR_uummdds, CR_hhmmsss, CR_energy, CR_ra, CR_dec = funcHM.read_TAascii(CRasciiFile, eThrPower, form='wiki')

#print CR_energy[:10]
#print CR_ra[:10]*180 / np.pi
#print CR_dec[:10]*180 / np.pi

#funcHM.drawSkymap(CR_ra, CR_dec)

# CR histogram
CR_hist_rad, CR_hist_binCenter, CR_binsizes = \
	funcHM.posDelta_hist(CR_energy, CR_ra, CR_dec)
#plt.plot(CR_hist_binCenter, CR_hist_rad, marker="o", label = 'CR events')
CR_hist_err = np.sqrt(CR_hist_rad)
plt.errorbar(CR_hist_binCenter, CR_hist_rad, xerr = CR_binsizes, yerr = CR_hist_err, ecolor='k', fmt = 'none', label = 'CR events')

# MC events
MC_uummdds, MC_hhmmsss, MC_energy, MC_ra, MC_dec = funcHM.read_TAascii(MCasciiFile, eThrPower, form='wiki')

#print MC_energy[:10]
#print MC_ra[:10]*180 / np.pi
#print MC_dec[:10]*180 / np.pi

# Seperate MC events to the number of CR events
MC_seperate = np.array(funcHM.splitArr([MC_energy, MC_ra, MC_dec], len(CR_energy)))

# Check first arrays set
#funcHM.drawSkymap(MC_seperate[0][1], MC_seperate[0][2])

# Draw histogram of position delta
total_hist = 0
total_hist_sq = 0
for i in range(len(MC_seperate)):
	hist_rad, hist_binCenter, binsizes = funcHM.posDelta_hist(MC_seperate[i][0], MC_seperate[i][1], MC_seperate[i][2])
	#print hist_rad
	#print hist_binCenter

	total_hist += hist_rad
	total_hist_sq += np.square(hist_rad)
hist_avg = total_hist / len(MC_seperate)
hist_error = np.sqrt(total_hist_sq / len(MC_seperate) - np.square(hist_avg))


#plt.plot(hist_binCenter, total_hist)
#plt.errorbar(hist_binCenter, hist_avg, hist_error, binsizes)
plt.plot(hist_binCenter, hist_avg, label = 'MC events', linestyle = ':')
plt.legend()
plot_title = str("E > %4.1f EeV" % eThrPower)
plt.title(plot_title)
plt.show()

#plt.plot(hist_binCenter, 1.*(CR_hist_rad - hist_avg) / hist_avg, label = 'MC events', linestyle = '-', color='black')
plt.errorbar(CR_hist_binCenter, 1.*(CR_hist_rad - hist_avg) / hist_avg, xerr = CR_binsizes, yerr = CR_hist_err / hist_avg, ecolor='k', fmt = 'none', label = 'CR events')
plt.title(plot_title)
plt.show()