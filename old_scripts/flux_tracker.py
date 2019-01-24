import numpy as np
import matplotlib.pyplot as plt
import csv
import operator
import time

#User inputs to find relevant csv file
system = str(raw_input('Planetary system: '))
planet = str(raw_input('Impacted planet: '))
ratio = float(raw_input('Ratio of kick to Keplerian velocity: '))
N_testparticle = int(raw_input('Number of test particles: '))
tp = int(raw_input('Test particle index: '))

if system == 'TRAPPIST-1':
	N_active = 8

#User specifies what waveband they'd like
lambda_1 = int(raw_input('Wavelength (lower bound) in Angstroms: '))
lambda_2 = int(raw_input('Wavelength (upper bound) in Angstroms: '))

#Constants
G = 6.67e-11

#Class for extracting parameters from file
class Satellite():
	def __init__ (self, data):
		self.data = data
		
	def __str__(self):
		return "Mass:\t\t\t\t" + str(self.getMass()) + "\nSemi Major Axis:\t\t" + str(self.getSemiMajorAxis()) + "\nEccentricity:\t\t\t" + str(self.getEccentricity()) + "\nInclination:\t\t\t" + str(self.getInclination()) + "\nMean Anomaly:\t\t\t" + str(self.getMeanAnomaly()) + "\nLongitude of Ascending Node:\t" + str(self.getLongAscNode()) + "\nArgument of Pericenter:\t\t" + str(self.getArgPericenter()) + "\nRadius:\t\t\t\t" + str(self.getRadius()) + "\n"
		
	def __lt__(self,other):
		if self.getSemiMajorAxis() < other.getSemiMajorAxis():
			return True
		if self.getSemiMajorAxis() == None:
			return True
		return False	
		
	def getDictionary(self):
		retDict = {}
		if self.getMass() != None:
			retDict['m'] = self.getMass()
		if self.getSemiMajorAxis() != None:
			retDict['a'] = self.getSemiMajorAxis()
		if self.getEccentricity() != None:
			retDict['e'] = self.getEccentricity()
		if self.getInclination() != None:
			retDict['inc'] = self.getInclination()
		if self.getMeanAnomaly() != None:
			retDict['M'] = self.getMeanAnomaly()
		if self.getLongAscNode() != None:
			retDict['Omega'] = self.getLongAscNode()
		if self.getArgPericenter() != None:
			retDict['omega'] = self.getArgPericenter()
		if self.getRadius() != None:
			retDict['r'] = self.getRadius()
		return retDict

	def getData(self):
		return self.data
	
	def getMass(self):
		if self.data[0]:
			return float(self.data[0])
		return None
	
	def getSemiMajorAxis(self):
		if self.data[1]:
			return float(self.data[1])
		return None
	
	def getEccentricity(self):
		if self.data[2]:
			return float(self.data[2])
		return None
	
	def getInclination(self):
		if self.data[3]:
			return float(self.data[3])
		return None
	
	def getMeanAnomaly(self):
		if self.data[4]:
			return float(self.data[4])
		return None
	
	def getLongAscNode(self):
		if self.data[5]:
			return float(self.data[5])
		return None
	
	def getArgPericenter(self):
		if self.data[6]:
			return float(self.data[6])
		return None
	
	def getRadius(self):
		if self.data[7]:
			return float(self.data[7])
		return None;

#Load data for appropriate planetary system
data = list(csv.reader(open('C:\Users\James\Documents\Warwick\Year 4\URSS\Planetary_systems\\'+str(system)+'_test.csv')))

satellites = {}
for i in data:
	if i[0] != "Object":
		satellites[i[0]] = Satellite(i[1:])

readKey=planet
if readKey in satellites:
	print satellites[readKey]
else:
	print readKey + " is not in satellites"

#Store parameters of star and planets in params
params = {}
for name, satellite in sorted(satellites.iteritems(), key=operator.itemgetter(1)):
		params.update({name:satellite.getDictionary()})

#Necessary parameters for later
a = params[planet]['a']
m_star = params['A']['m']
r_star = params['A']['r']

#Obtain time array for simulation
def orb_period(a,m_star):
	T = 2*np.pi*np.sqrt(np.power(a,3)/(G*m_star))
	return T

T = orb_period(a,m_star)

N_orbits = 1000
N_outputs = 1e6

times_secs = np.linspace(0,N_orbits*T,N_outputs) #For integration
times = [i / T for i in times_secs] #For plots

#Load the destinations file
merge_times = {}
with open('C:\Users\James\Documents\Warwick\Year 4\URSS\Output_files\\'+system+'\\'+planet+'\Merge_info_'+planet+'_'+str(ratio)+'_'+str(N_testparticle)+'.csv','r') as infile:
	reader = csv.reader(infile, delimiter = '\t')
	for row in reader:
		merge_times.update({row[0]:row[1]})

#Reduce the array to include just the active data
cut_off = float(merge_times[str(tp)])
print merge_times[str(17)]
print cut_off
times_cut = []
for i in times:
    if i < cut_off:
        times_cut.append(i)

times_cut = np.array(times_cut)

#Obtain array for distances of test particle from star during simulation
t0 = time.time()
print 'Extracting distances from file. Please wait...'
distances = np.genfromtxt('C:\Users\James\Documents\Warwick\Year 4\URSS\Output_files\\'+system+'\\'+planet+'\d_star_'+planet+'_'+str(ratio)+'_'+str(N_testparticle)+'.csv', skip_header=1, usecols=(tp-N_active), unpack=True, delimiter='\t', max_rows=len(times_cut))
print 'Done.'
t1 = time.time()
print t1-t0

#Check for zeros
mask = ~(distances == 0)
distances = distances[mask]
times_cut = times_cut[mask]

times_cut_secs = [i * T for i in times_cut]

#Obtain stellar spectrum
print 'Extracting stellar spectrum from file...'
spectrum = np.genfromtxt('C:\Users\James\Documents\Warwick\Year 4\URSS\Stellar_spectra\\'+system+'\lte26-5.0-0.0.AMES-dusty.7.dat.txt', delimiter = '\t', skip_header = 7)
print 'Done.'

#Reduce spectrum to desired waveband
wl = []
sfd = [] #Spectral flux density [W/m^2/m]
for i in spectrum:
	if i[0] >= lambda_1 and i[0] <= lambda_2:
		wl.append(i[0])
		sfd.append(i[1]) #This gives spectral flux density at surface of star

#Iterate through distances and calculate corresponding flux
irrad = [] #Irradiance [W/m^2]
print 'Normalising spectral flux density and calculating irradiance...'
for d in distances:
	sfd_norm = [i * np.square(r_star/d) for i in sfd] #Spectral flux density normalised for given distance
	irr = np.trapz(sfd_norm,wl)
	irrad.append(irr)
print 'Done'

irradiance = [i * 0.001 for i in irrad]

#Calculate the overall flux 
fluence = np.trapz(irradiance,times_cut_secs)
print fluence

plt.plot(times_cut,irradiance,'b.')
plt.xlabel('Time (orbits)')
plt.ylabel('Irradiance (J/m^2/s)')
plt.show()
