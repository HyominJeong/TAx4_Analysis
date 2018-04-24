import funcHM
import numpy as np

Energies = np.array([1, 1, 2, 3, 4, 5])
Ras = np.array([15, 25, 35, 45, 55, 55]) * np.pi / 180
Decs = np.array([85, 85, 75, 65, 55, 55]) * np.pi / 180

CR_events = [Energies, Ras, Decs]

result = funcHM.nOfEvts_RaDec(CR_events)

for i in range(len(result[0])):
	for j in range(len(result[0][i])):
		if result[0][i][j] > 0:
			print result[0][i][j], result[1][i][j]* 180 / np.pi
