import numpy as np

dis_vel_set = np.array([[15, 1100],
[97, 6700],
[32, 2400],
[145, 10700],
[50, 3100],
[122, 9900],
[58, 4300],
[91, 5300],
[120, 9000],
[93, 7500],
[158, 8900],
[64, 5300],
[145, 9600],
[61, 3300],
[103, 5100],
[46, 3600],
[34, 1800],
[185, 9500],
[20, 1200]])
print(dis_vel_set)
print(dis_vel_set[0])

for i in range(len(dis_vel_set)):
	print(dis_vel_set[i][0])

dis_vel_trans = dis_vel_set.transpose()
print dis_vel_trans
print dis_vel_trans[0][:10]