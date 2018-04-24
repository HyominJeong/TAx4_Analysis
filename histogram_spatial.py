import numpy as np
import matplotlib.pyplot as plt
import os, sys
import funcHM

def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False


infile = open(sys.argv[1])

events = []

thetas = []
eThresh = 1.*10**(1.) # energy threshold for analysis, in EeV
# Making events array
# events = [ [ E, [Ra,Dec], [X,Y,Z] ]
#            [          ...         ]
#            [                      ] ]
for line in infile:

	#print line.rstrip('\n').split()
	#print len(line.rstrip('\n').split())
	if len(line.rstrip('\n').split()) == 17:
		#head, date, hhmmss, xcore, ycore, s800, theta, phi, dummy1, dummy2, energy = line.rstrip('\n').split()
		head, yymmdd, hhmmss_usec, jday, lmst, energy, theta, dtheta, phi, dphi, ha, ra, dec, l, b, sgl, sgb = line.rstrip('\n').split()
		#print float(energy)
		#break
		if isfloat(energy):
			if float(energy)> eThresh:
				event = []
				event.append([float(energy)])
				event.append([float(ra),float(dec)])
				event.append(funcHM.radec2xyz(float(ra),float(dec)))
				events.append(event)
				#print event
				#break
				#energies.append(numpy.log10(float(energy)))
				#if numpy.log10(float(energy[:-1]))>18: and float(theta)>0 and float(phi)>0:
				#thetas.append(float(theta[:-1]))

# Get position deltas
print "Getting position deltas"
print "Total number of events is", len(events)
rad_deltas = []
for n in range(len(events)):
	for i in range(n):
		pos1 = events[i][2]
		pos2 = events[n][2]
		#print pos1, pos2, i, n
		len_delta = funcHM.posDelta(pos1, pos2)
		#print 2.*np.arcsin(0.5*len_delta)
		rad_deltas.append(2.*np.arcsin(0.5*len_delta))
		
#	if n == 100:
#		break


# Draw histogram of 
hist_rad_deltas, binEdge_rad_deltas = np.histogram(rad_deltas, bins=100)
bincenters = 0.5*(binEdge_rad_deltas[1:] + binEdge_rad_deltas[:-1])
binsizes = 0.5*(binEdge_rad_deltas[1:] - binEdge_rad_deltas[:-1])
menStd = np.sqrt(hist_rad_deltas)

plt.errorbar(bincenters, hist_rad_deltas, menStd, binsizes)
plt.show()


'''
# Draw energy distribution

bin_energy_log = [0.1*n for n in range(-1,25)]
hist_log_energy, binEdge_energy = numpy.histogram(energies, bins=bin_energy_log)

bincenters = 0.5*(binEdge_energy[1:] + binEdge_energy[:-1])
#print bincenters

menStd = numpy.sqrt(hist_log_energy)

#print len(bin_energy_log)
#print len(hist_log_energy[0])
#plt.plot(bincenters, hist_log_energy, drawstyle='steps-mid', color='black')
#binsizes = binEdge_energy[1:] - binEdge_energy[:-1]

#width = 0.05
#plt.scatter(bincenters, hist_log_energy, color='k')
plt.errorbar(bincenters, hist_log_energy, menStd, 0.05, color='k', ls='none', elinewidth=1)

plt.yscale("log")

plt.xlabel('log10(E/EeV)')
plt.ylabel('NumberOfEvents')

plt.show()

# Draw zenith angle distribution

bin_theta = [1.25*n for n in range(0,40)]

hist_theta, binEdge_theta = numpy.histogram(thetas, bin_theta)
print binEdge_theta
bincenters_theta = 0.5*(binEdge_theta[1:] + binEdge_theta[:-1])
menStd_theta = numpy.sqrt(hist_theta)

plt.errorbar(bincenters_theta, hist_theta, menStd_theta, 0.625, color='k', ls='none', elinewidth=1)

#plt.yscale("log")

plt.xlabel('Zenith angle(deg)')
plt.ylabel('NumberOfEvents')

plt.show()
'''

