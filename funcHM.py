from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from scipy.special import sph_harm
import numpy as np
import matplotlib.pyplot as plt
import os, sys

# Return position delta from Cartesian coordinate
# pos1, pos2 = [x,y,z] in numpy array
# output is float
# Hyomin Jeong
# 20180117
def posDelta(pos1,pos2):
	return(np.sqrt((pos1[0]-pos2[0])**2+(pos1[1]-pos2[1])**2+(pos1[2]-pos2[2])**2))

# Return [x,y,z] vector from Ra, Dec
# Ra, Dec in radians
# Output is [x,y,z] numpy array
# Hyomin Jeong
# 20180117
def radec2xyz(Ra, Dec):
	x = np.sin(0.5*np.pi-Dec)*np.cos(1.*Ra)
	y = np.sin(0.5*np.pi-Dec)*np.sin(1.*Ra)
	z = np.cos(0.5*np.pi-Dec)
	return(np.array([x,y,z]))

# Return [Ra, Dec] in radian from [X,Y,Z] float vector
# X, Y, Z can be numpy array
# Hyomin Jeong
# 20180130
def xyz2radec(X, Y, Z):
	R = np.sqrt(X**2 + Y**2 + Z**2)
	theta = np.arccos(Z / R)
	phi = np.arctan2(X, Y)
	Ra = phi
	Dec = 0.5 * np.pi - theta
	return Ra, Dec

# Return solid angle of each Ra, Dec bin
# ra_rad, dec_rad, delta_rad in radian unit
# Hyomin Jeong
# 20180130
def dohm(ra_rad, dec_rad, delta_rad):
	theta_rad = 0.5 * np.pi - dec_rad
	phi_rad = ra_rad

	# Canceled because '>' operator can not read numpy array
	#if theta_rad > np.pi:
	#	print "Strange Dec!", ra_rad, dec_rad, delta_rad

	result = delta_rad * \
			 (-np.cos(theta_rad+0.5*delta_rad) + np.cos(theta_rad-0.5*delta_rad))

 	# Canceled because '>' operator can not read numpy array
	#if result < 0:
	#	print "Strange Result!", ra_rad, dec_rad, delta_rad, result
	return result

# Return homogeneous Ra, Dec sets
# nOfEvts : integer
# Output is [[Ra sets], [Dec sets]] numpy array in radian unit
# Hyomin Jeong
# 20180131
def CRHradec(nOfRand, show=0):
	# Ra from 0 to 360
	ra_random = np.random.rand(nOfRand)*2*np.pi#*360.0*u.degree

	# Dec from -90 to 90
	#dec_random=(np.random.rand(nOfRand)*2 -1)*2*np.pi#*180.0-90.0)*u.degree
	theta_random=np.arccos(np.random.rand(nOfRand)*2-1)
	dec_random=(0.5 * np.pi - theta_random)# * 180 / np.pi *u.degree

	###########################################################
	# Show polar scatter plot of Ra, Dec = (0, 0) and (0, 90) #
	###########################################################
	if show == 1:
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

	return ra_random, dec_random

# Return position delta histogram from Ra, Dec sets
# Ra, Dec sets are array in radian unit
# Output is numpy histogram of position delta in radian unit
# Hyomin Jeong
# 20180131
def histPosDelta(Ras, Decs, nbin):
	print "Getting position deltas"
	print "Total number of events is", len(Ras)
	rad_deltas = []
	for n in range(len(Ras)):
		for i in range(n):
			pos1 = radec2xyz(float(Ras[i]),float(Decs[i]))
			pos2 = radec2xyz(float(Ras[n]),float(Decs[n]))
			#print pos1, pos2, i, n
			len_delta = posDelta(pos1, pos2)
			#print 2.*np.arcsin(0.5*len_delta)
			rad_deltas.append(2.*np.arcsin(0.5*len_delta))

	# Draw histogram of 
	hist_rad_deltas, binEdge_rad_deltas = np.histogram(rad_deltas, bins=nbin)
	return hist_rad_deltas, binEdge_rad_deltas

# Check input value is float or not
def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False


# Return energy, ra, dec sets from TA-anisotropy ascii
# Input: name of ascii file
# Output: numpy array of [ [Energy array], [Ra array], [Dec array] ]
# Energy in EeV unit
# Ra, Dec in radian unit
# Added collecting upper limit, 2018.03.15.
def read_TAascii(TAascii, eThr = 1, eLim = -1, form = 'anisotropy'):

	# Check energy limit
	if (eLim != -1) and (eThr > eLim):
		print "Invalid energy range!"
		sys.exit()

	print "Read TA ascii file", TAascii

	infile = open(TAascii)

	events = []

	thetas = []
	#eThresh = 1.*10**(eThrPower) # energy threshold for analysis, in EeV
	eThresh = eThr # energy threshold for analysis, in EeV
	#eLim = 1.*10**(eLimPower) # energy epper limit for analysis, in EeV
	
	# Make events array
	# result = [ [Energy array], [Ra array], [Dec array] ]

	yymmdds = np.array([])
	hhmmsss = np.array([])
	energies = np.array([])
	ras = np.array([])
	decs = np.array([])


	if form == 'wiki':
		print "################################################################################"
		print "# TA-wiki form ascii from                                                      #"
		print "# date, hhmmss, xcore, ycore, s800, theta_deg, phi_deg, dummy1, dummy2, energy #"
		print "################################################################################"
		for line in infile:

			#print line.rstrip('\n').split()
			#print len(line.rstrip('\n').split())
			if len(line.rstrip('\n').split()) == 10:
				yymmdd, hhmmss, xcore, ycore, s800, theta_deg, phi_deg, dummy1, dummy2, energy = line.rstrip('\n').split()
				#head, yymmdd, hhmmss_usec, jday, lmst, energy, theta, dtheta, phi, dphi, ha, ra, dec, l, b, sgl, sgb = line.rstrip('\n').split()
				if isfloat(energy):
					if float(energy)> eThresh:
						if (eLim == -1) or (float(energy) < eLim):
							#print yymmdd, hhmmss, theta, phi
							yymmdds = np.append(yymmdds, int(yymmdd[:-1]))
							hhmmsss = np.append(hhmmsss, float(hhmmss[:-1]))

							energies = np.append(energies, float(energy))
							# Convert theta, phi to ra, dec
							ha, dec = thetaphi2hadec(float(theta_deg[:-1])*np.pi/180, float(phi_deg[:-1])*np.pi/180)
							lmst = lst(yymmdd[:-1], hhmmss[:-1].split('.')[0])
							ra = lmst - ha
							#print yymmdd, hhmmss, theta[:-1], phi[:-1], ra, dec
							ras = np.append(ras, float(ra))
							decs = np.append(decs, float(dec))

	elif form == 'anisotropy':
		print "##################################################################################################"
		print "# TA-anisotropy form ascii from                                                                  #"
		print "# yymmdd, hhmmss_usec, jday, lmst, energy, theta, dtheta, phi, dphi, ha, ra, dec, l, b, sgl, sgb #"
		print "##################################################################################################"
		n = 0
		for line in infile:

			#print line.rstrip('\n').split()
			#print len(line.rstrip('\n').split())
			if len(line.rstrip('\n').split()) == 17:
				#date, hhmmss, xcore, ycore, s800, theta, phi, dummy1, dummy2, energy = line.rstrip('\n').split()
				head, yymmdd, hhmmss_usec, jday, lmst, energy, theta, dtheta, phi, dphi, ha, ra, dec, l, b, sgl, sgb = line.rstrip('\n').split()
				if isfloat(energy):
					# print energy, ra, dec
					if float(energy)> eThresh:

						yymmdds = np.append(yymmdds, yymmdd)
						hhmmsss = np.append(hhmmsss, hhmmss_usec)
						energies = np.append(energies, float(energy))
						ras = np.append(ras, float(ra))
						decs = np.append(decs, float(dec))

						'''
						# Check ra, dec is invalid
						# Confirmed from TA-anisotropy data
						# Hyomin Jeong
						# 20180205
						hhmmss = hhmmss_usec.split('.')[0]
						ha_, dec_ = thetaphi2hadec(float(theta), float(phi))
						lmst_ = lst(yymmdd, hhmmss)
						print ha, ha_, dec, dec_, lmst, lmst_
						if n >= 10:
							sys.exit()
						n += 1
						'''
					#else:
					#	print energy
			else:
				print "Check the input file"
				print "Analysis is stoped"
				sys.exit()

	else:
		print "Invalid input file."
		print "Analysis is stoped"
		sys.exit()

	print "Array of Energy, Ra, Dec is reported"
	if eLim != -1:
		print "Energy threshold is from",eThr, "to", eLim, "EeV"
	else:
		print "Energy threshold is ", eThr, "EeV"
	print "The number of events is", len(energies)

	return np.array(yymmdds), np.array(hhmmsss), np.array(energies), np.array(ras), np.array(decs)

# Draw CR events on Mollwiede skymap
# Input is numpy arrays of [Ra], [Dec] in radian units
def drawSkymap(Ra_rad, Dec_rad, dataName="No title"):

	c = SkyCoord(ra=Ra_rad * 180 / np.pi * u.degree, dec=Dec_rad * 180 / np.pi * u.degree, frame='icrs')

	plt.figure(figsize=(8,4.2))
	plt.subplot(111, projection="aitoff")
	
	if dataName == "No title":
		titleOfPlot = "No title"
	else:
		titleOfPlot = "Aitoff projection of " + dataName

	ra_rad = c.ra.wrap_at(180*u.deg).radian
	dec_rad = c.dec.radian

	plt.title(titleOfPlot, y = 1.08)
	plt.grid(True)
	#plt.plot(ra_rad, dec_rad, 'o', markersize=50, alpha=0.002)
	plt.plot(ra_rad, dec_rad, 'o', markersize=1, alpha=0.5)
	plt.subplots_adjust(top=0.95,bottom=0.0)
	plt.show()


# Return Local Sidereal Time(S.T. as below) of TA from day (yymmdd) and time (hhmmss)
# Local Sidereal Time in in radian unit
# CLF longitude and latitude is in degree unit
# Hyomin Jeong
# 20180202
def lst(yymmdd, hhmmss, CLF_lon_deg = -112.90875, CLF_lat_deg = 39.29693):
	CLF_lon_rad = CLF_lon_deg * np.pi / 180
	CLF_lat_rad = CLF_lat_deg * np.pi / 180

	yy = '{0:02d}'.format(int(yymmdd[-6:]) // 10000)
	MM = '{0:02d}'.format(int(yymmdd[-6:]) % 10000 // 100)
	dd = '{0:02d}'.format(int(yymmdd[-6:]) % 100)

	hh = '{0:02d}'.format(int(hhmmss) // 10000)
	mm = '{0:02d}'.format(int(hhmmss) % 10000 // 100)
	ss = '{0:02d}'.format(int(hhmmss) % 100)

	t_str = '20' + str(yy) + '-' + str(MM) + '-' + str(dd) + 'T' + str(hh) + ':' + str(mm) + ':' + str(ss)
	#print t_str

	t = Time(t_str, scale='utc', location = (CLF_lon_rad*u.degree, CLF_lat_rad*u.degree))
	
	#gsdt_rad = t.sidereal_time('apparent','greenwich').radian
	#lsdt_rad = gsdt_rad - CLF_lon_rad
	lsdt_rad = t.sidereal_time('mean', 'greenwich').radian + CLF_lon_rad
	
	while lsdt_rad < 0:
		lsdt_rad += 2 * np.pi
	while lsdt_rad > 2 * np.pi:
		lsdt_rad -= 2 * np.pi
	

	#print "Local Siderial Time is:", lsdt_rad
	return lsdt_rad


# Return [Ra, Dec] in radian unit from Azimuth and Zenith angle in radian unit
# Hyomin Jeong
# 20180202
def AzZe2radec(phi_rad, theta_rad, yymmdd, hhmmss, CLF_lon_deg = -112.90875, CLF_lat_deg = 39.29693):

	# Set CLF location
	CLF_lon_rad = CLF_lon_deg * np.pi / 180
	CLF_lat_rad = CLF_lat_deg * np.pi / 180

	El_rad = 0.5 * np.pi - theta_rad

	lsdt_rad = lst(yymmdd[-6:], hhmmss)

	# Convert phi to Azimuthal angle (TA use phi as countclockwise from East(X-axis))
	Az_rad = 0.5 * np.pi - phi_rad
	while Az_rad < 0:
		Az_rad += 2 * np.pi

	arr_lsdt = \
	np.array([ [ np.cos(lsdt_rad),   np.sin(lsdt_rad), 0 ],\
			   [ np.sin(lsdt_rad), - np.cos(lsdt_rad), 0 ],\
			   [                0,                  0, 1 ] ])
	
	arr_CLF_lat_rad = \
	np.array([ [   np.sin(CLF_lat_rad), 0, np.cos(CLF_lat_rad) ],\
			   [                     0, 1,                  0  ],\
			   [ - np.cos(CLF_lat_rad), 0, np.sin(CLF_lat_rad) ] ])

	arr_AzEl = \
	np.array([np.cos(El_rad) * np.cos(Az_rad),\
			  np.cos(El_rad) * np.sin(Az_rad),\
			  np.sin(El_rad)                 ])

	temp_arr = np.dot(arr_lsdt, arr_CLF_lat_rad)
	#print arr_lsdt, arr_CLF_lat_rad, temp_arr
	temp_arr = np.dot(temp_arr, arr_AzEl)
	#print temp_arr
	#print temp_arr

	ra = np.arctan2(temp_arr[1], temp_arr[0])
	dec = np.arcsin(temp_arr[2])
	return ra, dec



# Return hadec from phi theta
# Input: theta_rad = Zenith angle in radian unit
#        phi_rad = Azimuthe angle in radian unit, X-axis is east, countclockwise
# Output: ha, dec of celestial coordinate in radianunit
# Hyomin Jeong
# 20180206
def thetaphi2hadec(theta_rad, phi_rad, CLF_lon_deg = -112.90875, CLF_lat_deg = 39.29693):

	# Set CLF location to radian
	CLF_lon_rad = CLF_lon_deg * np.pi / 180
	CLF_lat_rad = CLF_lat_deg * np.pi / 180

	
	# Convert theta to elevation angle and phi to Azimuthal angle
	# (TA use phi as countclockwise from East(X-axis))
	El_rad = 0.5 * np.pi - theta_rad
	Az_rad = 0.5 * np.pi - phi_rad
	#lsdt_rad = lst(yymmdd[-6:], hhmmss)

	# atan2(-cos(azm),cos(lat)/tan(zen)-sin(lat)*sin(azm));
	ha = np.arctan2(-np.cos(phi_rad), \
					np.cos(CLF_lat_rad)/np.tan(theta_rad) - np.sin(CLF_lat_rad)*np.sin(phi_rad))

	while ha < 0:
		ha += 2 * np.pi

	dec = np.arcsin(np.cos(theta_rad)*np.sin(CLF_lat_rad)+np.sin(theta_rad)*np.cos(CLF_lat_rad)*np.sin(phi_rad));


	#ha = np.arctan2(temp_arr[1], temp_arr[0])
	#dec = np.arcsin(temp_arr[2])
	#print "Ha, Dec is :", ha, dec
	return ha, dec


# Split array(s) to random arrays consist of 'nOfComp' components
# Input: [ [Arr_1], [Arr_2], ..., [Arr_n] ]
#        Arrays have to have same number of components
# Output: [ [ [Arr_1[nOfComp]], [Arr_2[nOfComp]], ..., [Arr_n[nOfComp]] ],
#			  [Arr_1[nOfComp]], [Arr_2[nOfComp]], ..., [Arr_n[nOfComp]] ],
#										...								  ]
# Hyomin Jeong
# 20180207
def splitArr( Arrays, nOfComp):
	# Check form
	if len(Arrays[0]) == 1:
		print "This array cannot be separated"
		sys.exit()

	print len(Arrays), "array(s) will be separated to", len(Arrays[0]) / nOfComp, \
		  "arrays consist of", nOfComp, "Components"
	'''
	index = np.array(range(len(Arrays[0])))

	output_Arrays = []

	while len(index) >= nOfComp:
		print index
		pickIndex = np.random.choice(index, nOfComp, replace=False)
		#print pickIndex
		partOfArrays = []
		for n in range(len(Arrays)):
			
			Arr_pick = []
			for i in pickIndex:
				Arr_pick.append(Arrays[n][i])
			#print Arr_pick
			partOfArrays.append(Arr_pick)
		#print partOfArrays
		output_Arrays.append(partOfArrays)
		index = np.delete(index, pickIndex)
		print index, pickIndex
		print
	'''
	output_Arrays = []

	while len(Arrays[0]) >= nOfComp:
		index = np.array(range(len(Arrays[0])))
		pickIndex = np.random.choice(index, nOfComp, replace=False)
		#print pickIndex
		partOfArrays = []
		for n in range(len(Arrays)):
			
			Arr_pick = []
			for i in pickIndex:
				Arr_pick.append(Arrays[n][i])
			#print Arr_pick
			partOfArrays.append(Arr_pick)
			Arrays[n] = np.delete(Arrays[n], pickIndex)
		#print partOfArrays
		output_Arrays.append(partOfArrays)
		#index = np.delete(index, pickIndex)
		#print index, pickIndex
		#print
	return output_Arrays

# Combine two arrays consist of 'nOfComp' components
# Input: [ [Arr1_1], [Arr1_2], ..., [Arr1_n] ],
#		 [ [Arr2_1], [Arr2_2], ..., [Arr2_n] ],
#        Arrays have to have same number of sub rrays
# Output: [ [ [Arr_[nOfComp]], [Arr_2[nOfComp]], ..., [Arr_n[nOfComp]] ]
# Hyomin Jeong
# 20180226
def combArr( Array1, Array2):
	output_Arrays = []
	if len(Array1) != len(Array2):
		print "Wrong arrays"
	else:
		for i in range(len(Array1)):
			output_Array = np.append(Array1[i], Array2[i])
			output_Arrays.append(output_Array)

	return output_Arrays




# Get position deltas
# Input: [Energies], [Ras], [Decs] in numpy arrays or floats
# Output: position deltas in radian units
# Hyomin Jeong
# 20180207
def posDelta_hist(Energies, Ras, Decs, bin=np.arange(0., np.pi, np.pi/100)):
	print "Getting position deltas"
	print "Total number of events is", len(Energies)
	rad_deltas = []
	for n in range(len(Energies)):
		for i in range(n):
			pos1 = radec2xyz(Ras[i], Decs[i])
			pos2 = radec2xyz(Ras[n], Decs[n])
			#print pos1, pos2, i, n
			len_delta = posDelta(pos1, pos2)
			#print 2.*np.arcsin(0.5*len_delta)
			rad_deltas.append(2.*np.arcsin(0.5*len_delta))

	hist_rad_deltas, binEdge_rad_deltas = np.histogram(rad_deltas, bins=bin)
	bincenters = 0.5*(binEdge_rad_deltas[1:] + binEdge_rad_deltas[:-1])
	binsizes = 0.5*(binEdge_rad_deltas[1:] - binEdge_rad_deltas[:-1])
	menStd = np.sqrt(hist_rad_deltas)

	#plt.errorbar(bincenters, hist_rad_deltas, menStd, binsizes)
	#plt.show()
	return hist_rad_deltas, bincenters, binsizes

# Get zenith angle from latitude of detector
#						declination of direction
#						local siderial time
# Input: lat_det_rad in radian
#		 dec_rad in radian
#		 lst_rad in radian
def zenith_rad(lat_det_rad, dec_rad, lst_rad):
	el_rad = np.arcsin(\
						np.sin(lat_det_rad) * np.sin(dec_rad) +\
						np.cos(lat_det_rad) * np.cos(dec_rad) * np.cos(lst_rad))
	return 0.5 * np.pi - el_rad

# Get exposure rate from zenith angle and declination
def exposure_rate(zenith_arr, zenith_cut, dec_rad):
	exp = 0
	for zenith in zenith_arr:
		if zenith < zenith_cut:
			#exp += 1
			exp += np.cos(zenith)
	return 1.* exp / len(zenith_arr)

# Get differential exposure of declination
# Wrong, not used
def dif_exp(zenith_rad_arr, zenith_rad_cut, dec_rad):
	integ_exp = 0.
	
	for zenith_rad in zenith_rad_arr:
		#print zenith_rad, zenith_rad_cut
		if zenith_rad < zenith_rad_cut:
			#print zenith_rad
			integ_exp += np.sin(zenith_rad) * np.cos(zenith_rad)

	#print integ_exp			
	dif_exp = integ_exp / len(zenith_rad_arr) / np.cos(dec_rad)
	return dif_exp

# Find index from bin edge
def find_ind(bin_edge, value):
	for bin in bin_edge:
		if bin > value:
			return bin -1
			break

	
# Count CR events of each Ra and Dec bins
# Input: CR_evts = [[CR_Energy], [CR_Ra], [CR_Dec]]
#		 Ra_bin_size, Dec_bin_size = Bin sizes of Ra, Dec in radian unit
# Output: [[nOfEvts], ... ] = Number of events of each bin
#		  [[ra, dec], ... ] = Edge of each bin in radian units
# Ras and Decs are in radian unit
def nOfEvts_RaDec(CR_events, Ra_bin_size_deg = 10, Dec_bin_size_deg = 10):

	# Make Bins
	nOfRaBin  = int(360./( Ra_bin_size_deg))
	nOfDecBin = int(180./(Dec_bin_size_deg))

	print "Number of bin is :", nOfRaBin, " for Ra and", nOfDecBin, "for Dec =", nOfRaBin * nOfDecBin

	nOfEvts = np.zeros((nOfRaBin+1,nOfDecBin+1), dtype=np.int)
	#denOfEvts = np.zeros((nOfDecBin,nOfRaBin), dtype=np.float)

	# Fill Evts in each bin
	for i in range(len(CR_events[0])):
		#print ra_random[i], dec_random[i]
		Index_RaBin 	= int(CR_events[1][i] * 180 / np.pi // int(Ra_bin_size_deg))
		Index_DecBin 	= int((0.5 * np.pi - CR_events[2][i]) * 180 / np.pi // int(Dec_bin_size_deg))
		#print Index_RaBin, Index_DecBin
		nOfEvts[Index_RaBin,Index_DecBin] += 1
		#if Index_DecBin == 18 or Index_DecBin == -1:
		#	print CR_events[0][i], CR_events[1][i], CR_events[2][i]

	# Edge of each bin
	#bin_edge = np.zeros((nOfRaBin+1,nOfDecBin+1), dtype=np.dtype(object))
	#for i_ra in range(nOfRaBin):
	#	for j_dec in range(nOfDecBin):
	#		bin_edge[i_ra, j_dec] = \
	#			np.array([(i_ra)  * Ra_bin_size_deg  * np.pi / 180,\
	#					  (j_dec) * Dec_bin_size_deg * np.pi / 180 ])

	# Edge of each bin
	Ra_bin_edge = np.zeros((nOfRaBin+1,nOfDecBin+1), dtype=np.dtype(object))
	Dec_bin_edge = np.zeros((nOfRaBin+1,nOfDecBin+1), dtype=np.dtype(object))
	for i_ra in range(nOfRaBin):
		for j_dec in range(nOfDecBin):
			Ra_bin_edge[i_ra, j_dec] = \
				(i_ra)  * Ra_bin_size_deg  * np.pi / 180
			Dec_bin_edge[i_ra, j_dec] = \
				(j_dec) * Dec_bin_size_deg * np.pi / 180

	# Output
	return nOfEvts, Ra_bin_edge, Dec_bin_edge


# Make central position of bin edge array
# Example
# Input:	[ [0,0, 0.0], [1.0, 0.0], [2.0, 0.0],
#			  [0,0, 1.0], [1.0, 1.0], [2.0, 1.0],
#			  [0,0, 2.0], [1.0, 2.0], [2.0, 2.0] ]
# Output: 	[ [0.5, 0.5], [1.5, 0.5],
#			  [0.5, 1.5], [1.5, 1.5] ]
def bin_edge2cent(inputArray):
	result = []
	
	for i in range(len(inputArray)):
		print i
		print inputArray[i]
		bin_center = (inputArray[i][:-1] + inputArray[i][1:]) * 0.5
		result.append(bin_center)

	result_tr = []
	inputArray_tr = np.transpose(np.array(result))
	for i in range(len(inputArray_tr)):
		print i
		print inputArray_tr[i]
		bin_center = (inputArray_tr[i][:-1] + inputArray_tr[i][1:]) * 0.5
		result_tr.append(bin_center)
	
	return np.transpose(np.array(result_tr))

# Draw aitoff figure from given CR events
# On going
def drawFlux(CR_Ra_rad, CR_Dec_rad, binsize_deg = 20):
	###########################
	# Display # of Evts / Ohm #
	###########################

	# 0. Unify Ra form to (-pi to pi)
	CR_Ra_rad = RaConvt(CR_Ra_rad)

	# 1. Make Ohm bins
	delta_deg = binsize_deg # Size of bin
	delta_rad = delta_deg * np.pi / 180

	nOfRaBin = int(360./delta_deg)
	nOfDecBin = int(180./delta_deg)

	nOfEvts = np.zeros((nOfDecBin,nOfRaBin), dtype=np.int)
	denOfEvts = np.zeros((nOfDecBin,nOfRaBin), dtype=np.float)
	flux = np.zeros((nOfDecBin,nOfRaBin), dtype=np.float)

	# 2. Fill Evts in each bin
	for i in range(len(CR_Ra_rad)):
		#print ra_random[i], dec_random[i]
		Index_RaBin = int(((CR_Ra_rad[i] + np.pi) /np.pi * 180)%360 // int(delta_deg))
		Index_DecBin = int(((0.5 * np.pi - CR_Dec_rad[i])/np.pi * 180) // int(delta_deg))
		#print round(CR_Ra_rad[i]*180/np.pi), round(CR_Dec_rad[i]*180/np.pi), Index_RaBin, Index_DecBin, delta_deg
		nOfEvts[Index_DecBin,Index_RaBin] += 1

	#print nOfEvts
	#print np.sum(nOfEvts)

	# 3. Fill Ra and Dec of each bin
	#print nOfRaBin, nOfDecBin
	Pos_Ra = np.array([[(delta_deg * (n+0.5) - 180) * np.pi / 180 for n in range(nOfRaBin)] for i in range(nOfDecBin)])
	Edge_Ra = np.array([[(delta_deg * (n) - 180) * np.pi / 180 for n in range(nOfRaBin+1)] for i in range(nOfDecBin+1)])
	Pos_Dec = np.transpose([[(90. - (delta_deg * (n+0.5))) * np.pi / 180 for n in range(nOfDecBin)] for i in range(nOfRaBin)])
	Edge_Dec = np.transpose([[(90. - (delta_deg * (n))) * np.pi / 180 for n in range(nOfDecBin+1)] for i in range(nOfRaBin+1)])

	print "Edges of Ra\n", Edge_Ra[0] * 180 / np.pi
	print "Edgec of Dec\n", np.transpose(Edge_Dec)[0] * 180 / np.pi

	# 4. Calculate density and draw
	for i in range(len(nOfEvts)):
		#print funcHM.dohm((0.5*np.pi - Pos_Dec[i][0]), Pos_Ra[i][0], delta_rad), Pos_Dec[i][0], Pos_Ra[i][0]
		for j in range(len(nOfEvts[i])):
			#print Pos_Dec[i][j], Pos_Ra[i][j]
			denOfEvts[i,j] = 1.*nOfEvts[i,j] / (dohm(Pos_Ra[i][j], Pos_Dec[i][j], delta_rad))

			dir_exp = direc_exp(Pos_Dec[i][j], zen_cut_deg = 45)
			if dir_exp > 0:
				flux[i][j] = denOfEvts[i][j] / dir_exp			
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
	flux_disp = np.insert(flux, len(flux), 0, axis=0)
	flux_disp = np.insert(flux_disp, len(flux_disp[0]), 0, axis=1)

	#flux = np.zeros(denOfEvts_disp.shape)

	#print Edge_Ra.shape, Edge_Dec.shape, denOfEvts_disp.shape, Pos_Dec.shape

	CR_plt_alpha = 0.8

	fig1 = plt.figure()
	ax1 = fig1.add_subplot(111, projection="aitoff")
	ax1.grid(True)
	#p = ax1.pcolormesh(Pos_Ra, Pos_Dec, nOfEvts, alpha = 0.3)
	p = ax1.pcolormesh(Edge_Ra, Edge_Dec, denOfEvts_disp, alpha = 0.3)
	ax1.set_title("Density of events (# of evts / ohm)", y = 1.08)
	fig1.colorbar(p, orientation='horizontal')
	ax1.plot(CR_Ra_rad, CR_Dec_rad, '.', markersize=2, alpha=CR_plt_alpha, color='black')

	fig2 = plt.figure()
	ax2 = fig2.add_subplot(111,projection="aitoff")
	ax2.grid(True)
	p = ax2.pcolormesh(Edge_Ra, Edge_Dec, nOfEvts, alpha = 0.3)
	ax2.set_title("Number of events", y = 1.08)
	fig2.colorbar(p, orientation='horizontal')
	ax2.plot(CR_Ra_rad, CR_Dec_rad, '.', markersize=2, alpha=CR_plt_alpha, color='black')

	fig3 = plt.figure()
	ax3 = fig3.add_subplot(111, projection="aitoff")
	ax3.grid(True)
	#p = ax1.pcolormesh(Pos_Ra, Pos_Dec, nOfEvts, alpha = 0.3)
	p = ax3.pcolormesh(Edge_Ra, Edge_Dec, flux, alpha = 0.3)
	ax3.set_title("Differential flux of CR events (# of evts / ohm year km**2)", y = 1.08)
	fig3.colorbar(p, orientation='horizontal')
	ax3.plot(CR_Ra_rad, CR_Dec_rad, '.', markersize=2, alpha=CR_plt_alpha, color='black')
	'''
	# Print examples
	sample_dec_index = len(Pos_Ra) / 2
	print "N of Evt example"
	print nOfEvts[sample_dec_index],"\n", Pos_Ra[sample_dec_index] * 180 / np.pi, "\n", Pos_Dec[sample_dec_index] * 180 / np.pi
	print "Density :", denOfEvts_disp[sample_dec_index]
	print "DifFlux :", flux_disp[sample_dec_index]
	'''
	plt.show()

	return flux, Pos_Ra, Pos_Dec


# Return directional exposure from declination and zenith angle cut (optional)
def direc_exp(dec_rad, zen_cut_deg = 45):

	# Set parameters
	boarder_cut_efficiency = 0.8
	exposure_area = 700 			# km**2
	exposure_time = 7 				# years

	if zen_cut_deg == 45:
		direc_exp_file = open("/home/jhminie/Work/exposure/exp_dec_45.txt",'r')
	elif zen_cut_deg == 55:
		direc_exp_file = open("/home/jhminie/Work/exposure/exp_dec_55.txt",'r')
	else:
		sys.exit()

	#lines = direc_exp_file.readlines()
	
	dec_min = -25.
	dec_deg = dec_rad * 180 / np.pi

	exp_output = 0

	while dec_min < dec_deg:
		line = direc_exp_file.readline().rstrip()
		if len(line.split(' ')) == 2:
			#print line
			dec_min, exp_output = line.split(' ')
			dec_min = float(dec_min)
			#print dec_min
		else:
			break

	#print dec_min, exp_output, dec_rad, dec_deg

	direc_exp_file.close()



	return float(exp_output) * boarder_cut_efficiency * exposure_area * exposure_time

	
# Convert Ra range to -pi to pi in radian or -180 to 180 in degree or vice versa
# Input: numpy array of ra set (0 to 2 pi) or (0 to 360)
# Output: numpy array of ra set (-pi to pi) or (-180 to 180)
def RaConvt(RaSet, unit='rad', output='sym'):
	if output == 'sym':
		# 0 to 2 * pi -> -pi to pi
		RaSet = RaSet % (2 * np.pi)
		NofPi = RaSet // np.pi
		Remainder = RaSet % np.pi
		OutputRa = -1 * NofPi * np.pi + Remainder
	else:
		# -pi to pi -> 0 to pi
		RaSet = RaSet % (2 * np.pi)
	return OutputRa

# Calculate complex conjugate of Y_lm
# Input: [[flux], [Ra], [Dec], l, m]
# Output: [a_lm]
def a_lm(flux, Ra_rad, Dec_rad, l, m):
	Y_lm = sph_harm(m, l, Ra_rad, Dec_rad)
	Y_lm_conj = np.conj(Y_lm)

	# Calculate dohm to integrate
	delta_rad = Ra_rad[0][1] - Ra_rad[0][0]
	dohm_array = dohm(Ra_rad, Dec_rad, delta_rad)
	return np.sum(flux * Y_lm_conj * dohm_array)