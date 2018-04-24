import funcHM
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import matplotlib.pyplot as plt
import sys

CRasciiFile = sys.argv[1]
MCasciiFile = sys.argv[2]
eThrPower = float(sys.argv[3])

# CR events
CR_energy, CR_ra, CR_dec = funcHM.read_TAascii(CRasciiFile, eThrPower, form='wiki')

#print CR_energy[:10]
#print CR_ra[:10]*180 / np.pi
#print CR_dec[:10]*180 / np.pi

#funcHM.drawSkymap(CR_ra, CR_dec)

# CR histogram
CR_hist_rad, CR_hist_binCenter, CR_binsizes = \
	funcHM.posDelta_hist(CR_energy, CR_ra, CR_dec)
plt.plot(CR_hist_binCenter, CR_hist_rad, marker="o", label = 'CR events')

# MC events
MC_energy, MC_ra, MC_dec = funcHM.read_TAascii(MCasciiFile, eThrPower, form='wiki')

#print MC_energy[:10]
#print MC_ra[:10]*180 / np.pi
#print MC_dec[:10]*180 / np.pi

# Seperate MC events to the number of CR events for performance of computer
MC_seperate = np.array(funcHM.splitArr([MC_energy, MC_ra, MC_dec], len(CR_energy)))
print len(MC_seperate), "sets of MC events are read"
# Check first arrays set
#funcHM.drawSkymap(MC_seperate[0][1], MC_seperate[0][2])

# Draw histogram of position delta
total_hist = 0
total_hist_sq = 0

total_hist_HS = 0
total_hist_sq_HS = 0

for i in range(len(MC_seperate)):
	hist_rad, hist_binCenter, binsizes = funcHM.posDelta_hist(MC_seperate[i][0], MC_seperate[i][1], MC_seperate[i][2])


	# Add Hot spot
	Ohm_total = 2.10 	# Solid angle for 55 degree zenith cut
	Ohm_HS = 0.367 		# Solid angle for 20 degree hot spot (approximately)

	HS_ra_rad = 146.7 * np.pi / 180
	HS_dec_rad = 43.2 * np.pi / 180
	
	nOfEvt_HS = int(2.* len(MC_seperate[i][0]) / Ohm_total * Ohm_HS)
	HS_energy_arr = []
	HS_ra_arr = []
	HS_dec_arr = []

	print nOfEvt_HS, "Events are added for hot spot"
	#print len(MC_seperate[i])

	for j in range(nOfEvt_HS):
		HS_energy = 0.
		HS_ra = HS_ra_rad + 10 * np.pi / 180 * (np.random.rand() * 2 -1)
		HS_dec = HS_dec_rad + 10 * np.pi / 180 * (np.random.rand() * 2 -1)

		HS_energy_arr.append(HS_energy)
		HS_ra_arr.append(HS_ra)
		HS_dec_arr.append(HS_dec)

	# Append Hot spot events to MC events
	combMC = funcHM.combArr(MC_seperate[i], [HS_energy_arr, HS_ra_arr, HS_dec_arr])
		
	hist_rad_HS, hist_binCenter_HS, binsizes_HS = funcHM.posDelta_hist(combMC[0], combMC[1], combMC[2])
	#print hist_rad
	#print hist_binCenter

	total_hist += hist_rad
	total_hist_sq += np.square(hist_rad)

	total_hist_HS += hist_rad_HS
	total_hist_sq_HS += np.square(hist_rad_HS)	

hist_avg = total_hist / len(MC_seperate)
hist_error = np.sqrt(total_hist_sq / len(MC_seperate) - np.square(hist_avg))

hist_avg_HS = 1.* total_hist_HS / len(MC_seperate) * (1.* len(CR_energy) / (nOfEvt_HS + len(CR_energy)))**2
hist_error_HS = np.sqrt(total_hist_sq_HS / len(MC_seperate) - np.square(hist_avg_HS))


#plt.plot(hist_binCenter, total_hist)
#plt.errorbar(hist_binCenter, hist_avg, hist_error, binsizes)
plt.plot(hist_binCenter, hist_avg, label = 'MC events', linestyle = ':')
plt.plot(hist_binCenter_HS, hist_avg_HS, label = 'MC and HS events', linestyle = '-')
plt.legend()

plt.show()

plt.plot(hist_binCenter, CR_hist_rad - hist_avg, label = 'MC events', linestyle = '-', color='black')
plt.show()

