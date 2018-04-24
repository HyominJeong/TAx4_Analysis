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
infile2 = open(sys.argv[2])


energies = []
energies2 = []

thetas = []
thetas2 = []

# Data from input file 1
for line in infile:
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

# Data from input file 2
for line in infile2:
	if len(line.rstrip('\n').split()) == 10:
		date, hhmmss, xcore, ycore, s800, theta, phi, dummy1, dummy2, energy = line.rstrip('\n').split()
		#print float(energy)
		#break
		if isfloat(energy):
			if float(energy)>0: #and numpy.log10(float(energy[:-1]))<30:
				energies2.append(numpy.log10(float(energy)))
				#if numpy.log10(float(energy[:-1]))>18: and float(theta)>0 and float(phi)>0:
				thetas2.append(float(theta[:-1]))
		else:
			print date, event, energy, phi, theta

print numpy.average(energies), numpy.average(energies2)

# Draw energy distribution

bin_energy_log = [0.1*n for n in range(-1,25)]
hist_log_energy, binEdge_energy = numpy.histogram(energies, bins=bin_energy_log)
hist_log_energy2, binEdge_energy = numpy.histogram(energies2, bins=bin_energy_log)

bincenters = 0.5*(binEdge_energy[1:] + binEdge_energy[:-1])
#print bincenters

menStd = numpy.sqrt(hist_log_energy)
menStd2 = numpy.sqrt(hist_log_energy2)

#print len(bin_energy_log)
#print len(hist_log_energy[0])
#plt.plot(bincenters, hist_log_energy, drawstyle='steps-mid', color='black')
#binsizes = binEdge_energy[1:] - binEdge_energy[:-1]

#width = 0.05
#plt.scatter(bincenters, hist_log_energy, color='k')
plt.errorbar(bincenters, hist_log_energy, menStd, 0.05, color='k', ls='none', elinewidth=2)
plt.errorbar(bincenters, hist_log_energy2, menStd, 0.05, color='r', ls='none', elinewidth=1)

plt.yscale("log")

plt.xlabel('log10(E/EeV)')
plt.ylabel('NumberOfEvents')

plt.show()

# Draw zenith angle distribution

bin_theta = [1.25*n for n in range(0,40)]

hist_theta, binEdge_theta = numpy.histogram(thetas, bin_theta)
hist_theta2, binEdge_theta = numpy.histogram(thetas2, bin_theta)

#print binEdge_theta
bincenters_theta = 0.5*(binEdge_theta[1:] + binEdge_theta[:-1])
menStd_theta = numpy.sqrt(hist_theta)
menStd_theta2 = numpy.sqrt(hist_theta2)

plt.errorbar(bincenters_theta, hist_theta, menStd_theta, 0.625, color='k', ls='none', elinewidth=2)
plt.errorbar(bincenters_theta, hist_theta2, menStd_theta, 0.625, color='r', ls='none', elinewidth=1)

#plt.yscale("log")

plt.xlabel('Zenith angle(deg)')
plt.ylabel('NumberOfEvents')

plt.show()

