from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import matplotlib.pyplot as plt
import sys
import funcHM
#import matplotlib.patches as patches

nOfRand = 10000

# Old method to make Ra, Dec sets

# Ra from 0 to 360
ra_random = np.random.rand(nOfRand)*2*np.pi#*360.0*u.degree

# Dec from -90 to 90
#dec_random=(np.random.rand(nOfRand)*2 -1)*2*np.pi#*180.0-90.0)*u.degree
theta_random=np.arccos(np.random.rand(nOfRand)*2-1)
dec_random=(0.5 * np.pi - theta_random)# * 180 / np.pi *u.degree

###########################################################
# Show polar scatter plot of Ra, Dec = (0, 0) and (0, 90) #
###########################################################
fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, subplot_kw=dict(projection='polar'))

# Vernal equinox plot
theta = np.arctan2(np.sin(dec_random), np.cos(dec_random)*np.sin(ra_random))
r = np.sqrt((np.sin(ra_random)**2 + (np.sin(dec_random)**2*(np.cos(ra_random)**2))))

ax1.scatter(theta.flatten(), r.flatten(), s=0.5, alpha= 0.3, marker='.', color='black')
ax1.set_title("Center: Vernal equinox")

# Polar plot
theta=ra_random
r = np.cos(dec_random)

ax2.scatter(theta.flatten(), r.flatten(), s=0.5, alpha= 0.3, marker='.', color='black')
ax2.set_title("Center: North Pole")

# Y axis plot
theta = np.arctan2(np.sin(dec_random),- np.cos(dec_random)*np.sin(ra_random))
r = np.sqrt((np.sin(ra_random)**2 + (np.sin(dec_random)**2*(np.cos(ra_random)**2))))

ax3.scatter(theta.flatten(), r.flatten(), s=0.5, alpha= 0.3, marker='.', color='black')
ax3.set_title("Center: Y axis")


plt.show()

sys.exit()

'''
x_random = np.random.rand(nOfRand)*2 - 1
z_random = np.random.rand(nOfRand)*2 - 1
y_random = np.random.rand(nOfRand)*2 - 1

Ra, Dec = funcHM.xyz2radec(x_random, y_random, z_random)
ra_random = Ra * 180 / np.pi *u.degree
dec_random = Dec * 180 / np.pi *u.degree
'''
c = SkyCoord(ra=ra_random,dec=dec_random, frame='icrs')

#ra_deg = [0., 30., 60., 90., 120., 180.]*u.degree
#dec_deg = [0., 30., 0., -30., 90., -90.]*u.degree

#c = SkyCoord(ra=ra_deg,dec=dec_deg, frame='icrs')

ra_rad = c.ra.wrap_at(180*u.deg).radian
dec_rad = c.dec.radian

# Show Dec distribution
'''
plt.hist(dec_random/u.degree)
plt.show()
plt.hist(ra_random/u.degree)
plt.show()
'''



'''
plt.figure(figsize=(8,4.2))
plt.subplot(111, projection="aitoff")
plt.title("Aitoff projection of our random data")
plt.grid(True)
plt.plot(ra_rad, dec_rad, '.', markersize=2, alpha=0.3)
plt.subplots_adjust(top=0.95,bottom=0.0)
#plt.show()
'''
###########################
# Display # of Evts / Ohm #
###########################

# 1. Make Ohm bins
delta_deg = 20. # Size of bin
delta_rad = delta_deg * np.pi / 180

nOfRaBin = int(360./delta_deg)
nOfDecBin = int(180./delta_deg)

nOfEvts = np.zeros((nOfDecBin,nOfRaBin), dtype=np.int)
denOfEvts = np.zeros((nOfDecBin,nOfRaBin), dtype=np.float)

# 2. Fill Evts in each bin
for i in range(len(ra_random)):
	#print ra_random[i], dec_random[i]
	Index_RaBin = int((ra_random[i]/u.degree) // int(delta_deg))
	Index_DecBin = int((dec_random[i]/u.degree) // int(delta_deg))
	#print Index_RaBin, Index_DecBin
	nOfEvts[Index_DecBin,Index_RaBin] += 1

#print nOfEvts
#print np.sum(nOfEvts)

# 3. Fill Ra and Dec of each bin
#print nOfRaBin, nOfDecBin
Pos_Ra = [[(delta_deg * (n+0.5) - 180) * np.pi / 180 for n in range(nOfRaBin)] for i in range(nOfDecBin)]
Edge_Ra = [[(delta_deg * (n) - 180) * np.pi / 180 for n in range(nOfRaBin+1)] for i in range(nOfDecBin+1)]
Pos_Dec = np.transpose([[(90. - (delta_deg * (n+0.5))) * np.pi / 180 for n in range(nOfDecBin)] for i in range(nOfRaBin)])
Edge_Dec = np.transpose([[(90. - (delta_deg * (n))) * np.pi / 180 for n in range(nOfDecBin+1)] for i in range(nOfRaBin+1)])

print "Edges of Ra\n", Edge_Ra[0]
print "Edgec of Dec\n", np.transpose(Edge_Dec)[0]

# 4. Calculate density and draw
for i in range(len(nOfEvts)):
	#print funcHM.dohm((0.5*np.pi - Pos_Dec[i][0]), Pos_Ra[i][0], delta_rad), Pos_Dec[i][0], Pos_Ra[i][0]
	for j in range(len(nOfEvts[i])):
		#print Pos_Dec[i][j], Pos_Ra[i][j]
		denOfEvts[i,j] = 1.*nOfEvts[i,j] / (funcHM.dohm(Pos_Ra[i][j], Pos_Dec[i][j], delta_rad))
		#if (hm_func.dohm((0.5*np.pi - Pos_Dec[i][j]), Pos_Ra[i][j], delta_deg) < 0):
			#print Pos_Dec[i][j], Pos_Ra[i][j], hm_func.dohm((0.5*np.pi - Pos_Dec[i][j]), Pos_Ra[i][j], delta_deg)
			#break
		#print hm_func.dohm((0.5*np.pi - Pos_Dec[i][j]), Pos_Ra[i][0], delta_rad), Pos_Dec[i][j], Pos_Ra[i][j]
	#break

# Add garbage column and low at the end of value (End edge of Ra, Dec will not displayed)
denOfEvts_disp = np.insert(denOfEvts, len(denOfEvts), 0, axis=0)
denOfEvts_disp = np.insert(denOfEvts_disp, len(denOfEvts_disp[0]), 0, axis=1)
nOfEvts_disp = np.insert(nOfEvts, len(nOfEvts), 0, axis=0)
nOfEvts_disp = np.insert(nOfEvts_disp, len(nOfEvts_disp[0]), 0, axis=1)

print len(Edge_Ra), len(Edge_Dec), len(denOfEvts_disp)

fig1 = plt.figure()
ax1 = fig1.add_subplot(111, projection="aitoff")
ax1.grid(True)
#p = ax1.pcolormesh(Pos_Ra, Pos_Dec, nOfEvts, alpha = 0.3)
p = ax1.pcolormesh(Edge_Ra, Edge_Dec, denOfEvts_disp, alpha = 0.3)
ax1.set_title("Density of events (# of evts / ohm)", y = 1.08)
fig1.colorbar(p, orientation='horizontal')
ax1.plot(ra_rad, dec_rad, '.', markersize=1, alpha=1, color='black')

fig2 = plt.figure()
ax2 = fig2.add_subplot(111,projection="aitoff")
p = ax2.pcolormesh(Edge_Ra, Edge_Dec, nOfEvts, alpha = 0.3)
ax2.set_title("Number of events", y = 1.08)
ax2.grid(True)
fig2.colorbar(p, orientation='horizontal')
ax2.plot(ra_rad, dec_rad, '.', markersize=1, alpha=1, color='black')


plt.show()