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
times = [i / T for i in times_secs] #For plots and cut-offs

#Load the destinations file
merge_dict = {}
with open('C:\Users\James\Documents\Warwick\Year 4\URSS\Output_files\\'+system+'\\'+planet+'\Merge_info_'+planet+'_'+str(ratio)+'_'+str(N_testparticle)+'.csv','r') as infile:
    reader = csv.reader(infile, delimiter = '\t')
    for row in reader:
        merge_dict.update({row[0]:row[1]})

#Obtain stellar spectrum
spectrum = np.genfromtxt('C:\Users\James\Documents\Warwick\Year 4\URSS\Stellar_spectra\\'+system+'\lte26-5.0-0.0.AMES-dusty.7.dat.txt', delimiter = '\t', skip_header = 7)

#Reduce spectrum to desired waveband
wl = []
sfd = [] #Spectral flux density [ergs/cm^2/s/Angstrom]
for i in spectrum:
    if i[0] >= lambda_1 and i[0] <= lambda_2:
        wl.append(i[0])
        sfd.append(i[1]) #This gives spectral flux density at surface of star

#Create array of test particle indices to iterate through
tp_indices = np.arange(N_active,N_testparticle+N_active,1)

#Create fluence tracker file
with open('C:\Users\James\Documents\Warwick\Year 4\URSS\Output_files\\'+system+'\\'+planet+'\Fluences_'+planet+'_'+str(ratio)+'_'+str(N_testparticle)+'.csv', 'wb') as outfile:
        w = csv.writer(outfile, delimiter='\t')
        w.writerow(['Index', 'Fluence', 'Merge time'])

#Iterate through test particles
fluences = []
merge_times = []
for tp in tp_indices:
    
    #Reduce the array to include just the active data
    if str(tp) in merge_dict:
        cut_off = float(merge_dict[str(tp)])
        
        times_cut = []
        for i in times:
            if i < cut_off:
                times_cut.append(i)
    else:
        times_cut = times
        cut_off = 1000.
    
    times_cut = np.array(times_cut)
    
    t0 = time.time()
    #Obtain distances array from file
    distances = np.genfromtxt('C:\Users\James\Documents\Warwick\Year 4\URSS\Output_files\\'+system+'\\'+planet+'\d_star_'+planet+'_'+str(ratio)+'_'+str(N_testparticle)+'.csv', skip_header=1, usecols=(tp-N_active), unpack=True, delimiter='\t', max_rows=len(times_cut))
    
    #Check for zeros
    mask = ~(distances == 0)
    distances = distances[mask]
    times_cut = times_cut[mask]
    
    times_cut_secs = [i * T for i in times_cut]
    
    #Iterate through distances and calculate corresponding flux
    irradiance = [] #Irradiance [ergs/cm^2/s]
    for d in distances:
        sfd_norm = [i * np.square(r_star/d) for i in sfd] #Spec flux density normalised for given distance
        irr = np.trapz(sfd_norm,wl)
        irradiance.append(irr * 0.001) #Converts from ergs/cm^2/s to J/m^2/s
    
    #Calculate the fluence [J/m^2]
    fluence = np.trapz(irradiance,times_cut_secs)
    
    print 'Test particle: ' + str(tp) + ', Fluence: ' + str(fluence)
    
    fluences.append(fluence)
    merge_times.append(cut_off)
    
    with open('C:\Users\James\Documents\Warwick\Year 4\URSS\Output_files\\'+system+'\\'+planet+'\Fluences_'+planet+'_'+str(ratio)+'_'+str(N_testparticle)+'.csv', 'a') as outfile:
        w = csv.writer(outfile, delimiter='\t')
        w.writerow([tp, fluence, cut_off])
    
    t1 = time.time()
    
    print t1 - t0

#Create array for phi, angle of ejection
phi = np.linspace(0.,2.*np.pi,N_testparticle)

#Plot fluences vs time taken as function of phi
cm = plt.cm.get_cmap('RdYlBu')
sc = plt.scatter(merge_times, fluences, c = phi, s = 35, cmap = cm)
plt.colorbar(sc)
plt.savefig('C:\Users\James\Documents\Warwick\Year 4\URSS\Figures\Flux\Fluences_'+planet+'_'+str(ratio)+'_'+str(N_testparticle)+'.png')

plt.show()