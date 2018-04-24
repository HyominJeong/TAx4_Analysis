import funcHM
import numpy as np

inputArr = [ range(15), range(15, 30), range(30, 45)]
nOfComp = 4

print inputArr

output_Arrs = funcHM.splitArr(inputArr, nOfComp)

print np.array(output_Arrs)