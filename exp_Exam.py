import funcHM
import matplotlib.pyplot as plt
import numpy as np

lst_arr = np.arange(0, 2 * np.pi, np.pi / 1000)
#print lst_arr

# Set properties
CLF_lat_deg = 39.29693
Zen_angle_cut_deg = 55.0

CLF_lat_rad = CLF_lat_deg * np.pi / 180
Zen_angle_cut_rad = Zen_angle_cut_deg * np.pi / 180

'''
# Draw zenith angle of each declination
dec_deg_arr = [90.0, 70.0, 50.0, 30.0, 10.0, -10.0]
Zen_angle_cut_deg_arr = [Zen_angle_cut_deg for n in range(len(lst_arr))]

for dec_deg in dec_deg_arr:
	dec_rad = dec_deg * np.pi / 180
	theta_deg = []
	for lst in lst_arr:
		theta_deg.append(180. / np.pi * funcHM.zenith_rad(CLF_lat_rad, dec_rad, lst))
	txt_label = "declination : " + str(dec_deg)

	#exp_arr.append(funcHM.exposure_rate(theta_deg, Zen_angle_cut_deg))
	#print exp_arr[-1]

	plt.plot(lst_arr, theta_deg, label = txt_label)

plt.plot(lst_arr, Zen_angle_cut_deg_arr, 'k:', label = "Zenith angle cut" + str(Zen_angle_cut_deg) + "deg")

plt.legend()
plt.xlabel("Local sidereal time (radian)")
plt.ylabel("Zenith angle (deg)")
plt.show()
'''

# Draw exposure rate
dec_deg_arr = np.arange(-25.0, 90.0, 0.1)
exp_output = open('exp_dec_55.txt', 'w')
exp_arr = []
for dec_deg in dec_deg_arr:
	dec_rad = dec_deg * np.pi / 180
	theta_rad = []
	for lst in lst_arr:
		theta_rad.append(funcHM.zenith_rad(CLF_lat_rad, dec_rad, lst))
	#exp_arr.append(funcHM.exposure_rate(theta_rad, Zen_angle_cut_rad))
	exp_dec = funcHM.exposure_rate(theta_rad, Zen_angle_cut_rad, dec_rad)
	exp_arr.append(exp_dec)

	print>>exp_output, dec_deg, exp_dec
exp_output.close()
exp_arr = np.array(exp_arr)

# Draw exposure rate for 45 deg zenith cut

# Set properties
CLF_lat_deg = 39.29693
Zen_angle_cut_deg = 45.0

CLF_lat_rad = CLF_lat_deg * np.pi / 180
Zen_angle_cut_rad = Zen_angle_cut_deg * np.pi / 180
dec_deg_arr = np.arange(-25.0, 90.0, 0.1)
exp_output = open('exp_dec_45.txt', 'w')
exp_arr_45 = []
for dec_deg in dec_deg_arr:
	dec_rad = dec_deg * np.pi / 180
	theta_rad = []
	for lst in lst_arr:
		theta_rad.append(funcHM.zenith_rad(CLF_lat_rad, dec_rad, lst))
	#exp_arr.append(funcHM.exposure_rate(theta_rad, Zen_angle_cut_rad))
	exp_dec = funcHM.exposure_rate(theta_rad, Zen_angle_cut_rad, dec_rad)
	exp_arr_45.append(exp_dec)

	print>>exp_output, dec_deg, exp_dec
exp_output.close()
exp_arr_45 = np.array(exp_arr_45)

#print dec_deg_arr
#print exp_arr
plt.plot(dec_deg_arr, exp_arr * 700 * 7 * 0.8, 'g--')
plt.plot(dec_deg_arr, exp_arr_45 * 700 * 7 * 0.8, 'g-')

plt.show()
