import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('C:\Users\James\Documents\Warwick\Year 4\URSS\Stellar_spectra\TRAPPIST-1\lte26-5.0-0.0.AMES-dusty.7.dat.txt', delimiter = '\t', skip_header = 7)

wl = []
flux = []
for i in data:
	wl.append(i[0])
	flux.append(i[1])

plt.figure(figsize=(0.4,0.4))

plt.plot(wl,flux)
plt.xlabel('Wavelength (Angstroms)')
plt.ylabel('Flux (ergs/cm^2/s/Angstrom)')
plt.show()

'''
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt

data = ascii.read('/home/james/Documents/URSS/Stellar spectra/lte26-5.0-0.0.AMES-dusty.7')

wl = []
for string in data['col1']:
	value = string.replace('D','E')
	value = float(value)
	wl.append(value)

flux = []
for string in data['col2']:
	value = string.replace('D','E')
	value = 10*float(value) #convert to units of ergs/s/cm^2/A
	flux.append(value)

plt.plot(wl,flux,'b.')

plt.show()
'''
'''
import pyfits
import numpy as np
import matplotlib.pyplot as plt

dat = pyfits.open('/home/james/Documents/URSS/Stellar spectra/ckp00_3500.fits')	

tbdata = dat[1].data

wl = []
flux = []
for i in tbdata:
	wl.append(i[0])
	flux.append(i[11]*3.336e-19*np.square(i[0])*np.power(4*np.pi,-1))

plt.plot(wl,flux)
plt.show()
'''


