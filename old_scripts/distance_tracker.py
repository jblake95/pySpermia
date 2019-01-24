import csv
import matplotlib.pyplot as plt
import numpy as np

system = str(raw_input('Planetary system: '))
planet = str(raw_input('Impacted planet: '))
ratio = float(raw_input('Ratio of kick to Keplerian velocity: '))
N_testparticle = int(raw_input('Number of test particles: '))
N_active = int(raw_input('Number of active bodies: '))
N_outputs = int(raw_input('Number of outputs: '))
#N_orbits = raw_input('Number of orbits: ')
tp = int(raw_input('Test particle index: '))
#key = str(raw_input('Target planet: '))

#Orbital period
def orb_period(a):
	T = 2*np.pi*np.sqrt(np.power(a,3)/(G*m_star))
	return T

#Calculate orbital period of innermost planet for integration
a_innermost = 1.67e9
m_star = 1.6e29
G = 6.67e-11
T = orb_period(a_innermost)

times_full = np.linspace(0,1000.*T,N_outputs)

distances = np.genfromtxt('/home/james/Documents/URSS/Output_files/'+system+'/d_star_'+planet+'_'+str(ratio)+'_'+str(N_testparticle)+'.csv', skip_header=1, usecols=(tp-N_active), unpack=True, delimiter='\t')

times = [i / T for i in times_full]
distances = [i / (1.5e11) for i in distances]

print len(times)
print len(distances)

plt.plot(times,distances,'b')

plt.xlabel('Time (orbits of planet b)')
plt.ylabel('Relative distance (AU)')

plt.show()
