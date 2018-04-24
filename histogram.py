import numpy
import matplotlib.pyplot as plt
import os, sys

def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False


infile = open(sys.argv[1])

energies = []

thetas = []


for line in infile:

	#print line.rstrip('\n').split()
	#print len(line.rstrip('\n').split())
	if len(line.rstrip('\n').split()) == 11:
		head, date, hhmmss, xcore, ycore, s800, theta, phi, dummy1, dummy2, energy = line.rstrip('\n').split()
		#print float(energy)
		#break
		if isfloat(energy):
			if float(energy)>0: #and numpy.log10(float(energy[:-1]))<30:
				energies.append(numpy.log10(float(energy)))
				#if numpy.log10(float(energy[:-1]))>18: and float(theta)>0 and float(phi)>0:
				thetas.append(float(theta[:-1]))
		else:
			print date, event, energy, phi, theta

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

