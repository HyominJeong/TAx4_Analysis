import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy import units

#filename = 'WISE_W4SNRge3_and_W4MPRO_lt_6.0_RADecl_nohdr.dat'
#data = ascii.read(filename)
#coords = SkyCoord(ra=data['ra'], dec=data['dec'], unit='degree')

# Make CR events
nOfRand = 100000

# Ra from 0 to 360
ra_random = np.random.rand(nOfRand)*360.0*units.degree

# Dec from -90 to 90
#theta_random=np.arccos(np.random.rand(nOfRand)*2-1)
#dec_random=(0.5 * np.pi - theta_random) * 180 / np.pi *units.degree
dec_random = (np.random.rand(nOfRand)*2-1)*90.0*units.degree

coords = SkyCoord(ra=ra_random,dec=dec_random, frame='icrs')

ra = coords.ra.wrap_at(180 * units.deg).radian
dec = coords.dec.radian


color_map = plt.cm.Spectral_r
fig = plt.figure(figsize=(6, 4))
fig.add_subplot(111, projection='aitoff')
image = plt.hexbin(ra, dec, cmap=color_map,
                   gridsize=45, mincnt=1, bins='log')

plt.xlabel('R.A.')
plt.ylabel('Decl.')
plt.grid(True)
plt.colorbar(image, spacing='uniform', extend='max')
plt.show()