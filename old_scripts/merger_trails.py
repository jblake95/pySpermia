import numpy as np
import matplotlib.pyplot as plt
import csv
import operator

system = str(raw_input('Planetary system: '))
planet = str(raw_input('Impacted planet: '))
ratio = float(raw_input('Ratio of kick to Keplerian velocity: '))
N_testparticle = int(raw_input('Number of test particles: '))

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

#Obtain time array for simulation
def orb_period(a,m_star):
    T = 2*np.pi*np.sqrt(np.power(a,3)/(G*m_star))
    return T

T = orb_period(a,m_star)

N_orbits = 1000
N_outputs = 1e6

times = np.linspace(0,N_orbits*T,N_outputs)
t_orbits = [i / T for i in times] 

merge_times = []
planets = []
with open('C:\Users\James\Documents\Warwick\Year 4\URSS\Output_files\\'+system+'\\'+planet+'\Merge_info_'+planet+'_'+str(ratio)+'_'+str(N_testparticle)+'.csv', 'rb') as infile:
    reader = csv.reader(infile, delimiter='\t')
    for row in reader:
        merge_times.append(row[1])
        planets.append(row[2])

merge_times = merge_times[1:]
planets = planets[1:]

N_testparticle = float(N_testparticle)

count = 0
[n_A,n_b,n_c,n_d,n_e,n_f,n_g,n_h,n_X,n_total] = [0,0,0,0,0,0,0,0,0,0]
t_outputs = []
[N_A,N_b,N_c,N_d,N_e,N_f,N_g,N_h,N_X,N_total] = [[],[],[],[],[],[],[],[],[],[]]

while count < len(t_orbits):
    
    for i in range(len(merge_times)):
        
        if float(merge_times[i]) < t_orbits[count]:
            if planets[i] == 'A':
                n_A += 1
            elif planets[i] == 'b':
                n_b += 1
            elif planets[i] == 'c':
                n_c += 1
            elif planets[i] == 'd':
                n_d += 1
            elif planets[i] == 'e':
                n_e += 1
            elif planets[i] == 'f':
                n_f += 1
            elif planets[i] == 'g':
                n_g += 1
            elif planets[i] == 'h':
                n_h += 1
            elif planets[i] == 'X':
                n_X += 1
    
    n_total = n_A + n_b + n_c + n_d + n_e + n_f + n_g + n_h + n_X
    
    N_A.append((n_A / N_testparticle) * 100)
    N_b.append((n_b / N_testparticle) * 100)
    N_c.append((n_c / N_testparticle) * 100)
    N_d.append((n_d / N_testparticle) * 100)
    N_e.append((n_e / N_testparticle) * 100)
    N_f.append((n_f / N_testparticle) * 100)
    N_g.append((n_g / N_testparticle) * 100)
    N_h.append((n_h / N_testparticle) * 100)
    N_X.append((n_X / N_testparticle) * 100)
    N_total.append((n_total / N_testparticle) * 100)
    
    print count, N_total[count]
    
    [n_A,n_b,n_c,n_d,n_e,n_f,n_g,n_h,n_X,n_total] = [0,0,0,0,0,0,0,0,0,0]
    
    t_outputs.append(t_orbits[count])
    
    count += 1

plt.plot(t_outputs,N_A,'orange')
plt.plot(t_outputs,N_b,'firebrick')
plt.plot(t_outputs,N_c,'magenta')
plt.plot(t_outputs,N_d,'b')
plt.plot(t_outputs,N_e,'r')
plt.plot(t_outputs,N_f,'g')
plt.plot(t_outputs,N_g,'c')
plt.plot(t_outputs,N_h,'purple')
plt.plot(t_outputs,N_X,'k--')
plt.plot(t_outputs,N_total,'k:')

plt.xscale('log')
plt.yscale('log')

plt.xlim(0.09,1010)

plt.xlabel('Time (orbits of planet g)')
plt.ylabel('Cumulative N (%)')

plt.show()