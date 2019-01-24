import rebound
import numpy as np
import matplotlib.pyplot as plt
import math
import time
import csv
import operator
from collections import Counter

#Constants
m_sun = 1.989e30
m_jupiter = 1.898e27
m_earth = 5.972e24
a_sun_earth = 1.5e11
deg_to_rad = np.pi/180.
r_sun = 6.957e8
r_earth = 6.371e6
G = 6.674e-11

#User input
system = str(raw_input('Planetary system: '))
planet = str(raw_input('Impacted planet: '))
ratio = float(raw_input('Ratio of kick to Keplerian velocity: '))
N_testparticle = int(raw_input('Number of test particles: '))

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
data = list(csv.reader(open('/home/james/Documents/URSS/Planetary_systems/'+str(system)+'_test.csv')))

satellites = {}
for i in data:
	if i[0] != "Object":
		satellites[i[0]] = Satellite(i[1:])

readKey=planet
if readKey in satellites:
	print satellites[readKey]
else:
	print readKey + " is not in satellites"

#Require mean anomaly and argument of pericenter of chosen planet
[M_planet,omega_planet] = [satellites[planet].getMeanAnomaly(),satellites[planet].getArgPericenter()]

#For circular and coplanar orbits, need the following well-defined functions
def MeanLongitude(Omega,omega,M):
	l = Omega + omega + M
	return l

def LongPericenter(Omega,omega):
	pomega = Omega + omega
	return pomega

#Set up a rebound simulation
params = {}
def setupSimulation():
	sim = rebound.Simulation()
	sim.G = 6.674e-11
	#Iterate through the bodies
	for name, satellite in sorted(satellites.iteritems(), key=operator.itemgetter(1)):
		p = satellite.getDictionary() #Obtain dictionary of parameters passed via csv file
		
		#Check that input angles are in range [-pi,pi]
		if 'omega' in p:
			if p['omega'] > np.pi:
				p['omega'] -= 2*np.pi
		if 'M' in p:
			if p['M'] > np.pi:
				p['M'] -= 2*np.pi
		
		#Without generalisation, need mean anomaly and argument of pericenter of chosen planet to be zero
		if name == planet:
			[p['M'],p['omega']] = [0.,0.] #Set mean anomaly and argument of pericenter of chosen planet to zero
		elif name == 'A':
			None
		else:
			#p['M'] = np.random.uniform(-np.pi,np.pi) #Other simulations find that the starting mean anomalies have little effect, so randomise
			p['M'] -= M_planet #Maintain original differences between mean anomalies - randomise later
			p['omega'] -= omega_planet #Maintain original differences between arguments of periastron
			
			#Ensure updated values for mean anomaly and argument of pericenter lie within range [-pi,pi]
			if p['omega'] > np.pi:
				p['omega'] -= 2*np.pi 
			elif p['omega'] < -np.pi:
				p['omega'] += 2*np.pi
			if p['M'] > np.pi:
				p['M'] -= 2*np.pi 
			elif p['M'] < -np.pi:
				p['M'] += 2*np.pi
		
		#Find mean longitude and longitude of pericenter for each planet
		if 'omega' in p:
			l = MeanLongitude(p['Omega'],p['omega'],p['M'])
			pomega = LongPericenter(p['Omega'],p['omega'])
		
		    #Finalise dictionary of parameters before adding to rebound simulation
			del p['omega']
			del p['M']
			p.update({'l':l}) #Well defined for circular orbits
			p.update({'pomega':pomega}) #Well defined for coplanar orbits
		
		#Add the body to the rebound simulation
		params.update({name:p}) #Create dictionary of parameters needed later
		sim.add(**p) #Add the star/planets to the rebound simulation
	return sim

sim = setupSimulation()

#Set the N_active variable of rebound to the number of active (massive) particles in the simulation
sim.N_active = sim.N

#Store parameters of star and chosen planet for later
[m_star,r_star] = [params['A']['m'],params['A']['r']]
[m,a,e,inc,f,Omega,omega,r] = [params[planet]['m'],params[planet]['a'],params[planet]['e'],params[planet]['inc'],params[planet]['l'],params[planet]['Omega'],params[planet]['pomega'],params[planet]['r']]

#Calculate the definitions of a close encounter for the star and planets
def Hill_radius(a,m):
	r_H = a*np.power(m/(3*m_star),1./3.) #Hill radius function
	return r_H

r_Hill = {} #To contain the Hill radius definitions for close encounters
r_merge = {}
pdict = {} #To set up a dictionary which will be used to keep track of the distances throughout integration
for key in sorted(params):
	if key == 'A':
		r_H = 2*r_star
		r_m = 2*r_star
	else:
		[a_p,m_p,r_p] = [params[key]['a'],params[key]['m'],params[key]['r']]
		r_H = Hill_radius(a_p,m_p)
		r_m = 2*r_p
	r_Hill.update({key:r_H})
	r_merge.update({key:r_m})
	pdict.update({key:{}})

tp_indices = np.arange(sim.N_active,N_testparticle+sim.N_active,1)

for key in pdict:
    for i in tp_indices:
	    pdict[key].update({i:[]})

#Plot the orbits of the planetary system
#fig = rebound.OrbitPlot(sim, slices = False, color = False, periastron = False, trails = False)
#fig = rebound.OrbitPlot(sim, slices = True, color = True, periastron = True, limz = 1e10)
#fig.savefig('/home/james/Documents/URSS/Figures/TRAPPIST-1_orbits.png')
#plt.show()
#time.sleep(2000)

#Calculate the new parameters for the ejecta orbits
def a_ratio(phi):
	a_ratio = 1-2*ratio*np.sin(theta)*np.sin(phi)-np.square(ratio)
	return a_ratio

def post_a(phi):
	a_rat = a_ratio(phi)
	a_post = a / a_rat
	return a_post - 3*r_Hill[planet]

def ang_mom_ratio(phi):
	H = 1+2*ratio*np.sin(theta)*np.sin(phi)+np.square(ratio)*(np.square(np.cos(theta))+np.square(np.sin(theta)*np.sin(phi)))
	return H

def post_e(phi):
	a_post = a / a_ratio(phi)
	H = ang_mom_ratio(phi)
	e_post_squared = 1. - H*(a / a_post)
	e_post = np.sqrt(e_post_squared)
	return e_post

def post_f(phi):
    num = np.absolute(1. + ratio*np.sin(phi))*(1./np.tan(phi))
    den = 2. + ratio*np.sin(phi)
    tan_f_post = num / den
    f_post = np.arctan(tan_f_post)
    if phi > np.pi:
	    f_post -= np.pi
    if f_post < -np.pi:
	    f_post += 2*np.pi #Ensures that all f_post lie in range [-pi,pi]
    return f_post

def post_omega(phi):
	f_post = post_f(phi)
	omega_post = -f_post 
	return omega_post

def post_Omega():
	Omega_post = f
	return Omega_post

def post_pomega(phi):
	Omega_post = post_Omega()
	omega_post = post_omega(phi)
	pomega_post = Omega_post + omega_post + 0.1
	return pomega_post

#Obtain arrays of parameters for range of phi
theta = 90.*np.pi / 180. #Coplanar orbits
phi = np.linspace(0.,2.*np.pi,N_testparticle)

'''
#Test individual orbits
a_p=post_a(3)
e_p=post_e(3)
f_p=post_f(3)
omega_p=post_omega(3)
sim.add(a=a_p,e=e_p,inc=0.,M=f_p,Omega=0,omega=omega_p)
fig = rebound.OrbitPlot(sim, slices = True, color = True, periastron = True, limz = 1e10)
fig.show()
time.sleep(30)
'''

#Add the test particles to the rebound simulation
for i in phi:
	a_p = post_a(i)
	e_p = post_e(i)
	f_p = post_f(i)
	pomega_p = post_pomega(i)
	sim.add(a=a_p,e=e_p,inc=0.,f=f_p,Omega=0,pomega=pomega_p)

#Orbital period required for integration
def orb_period(a):
	T = 2*np.pi*np.sqrt(np.power(a,3)/(G*m_star))
	return T

#fig = rebound.OrbitPlot(sim, slices = True, color = True, periastron = False, limz = 1e10)
#fig.show()
#time.sleep(2000)

'''
########################################################################
###### WHFAST integrator used to investigate test paticle orbits #######
########################################################################

#Switch to centre of mass frame
sim.move_to_com()
sim.integrator = "whfast"
	
#Determine suitable timestep
a_innermost = params['b']['a']
T = orb_period(a_innermost)
sim.dt = T / 100

#Run the simulation for 200 orbits, keeping store of the positions of the particles 10 times during the interval
N_orbit = 100
t_max = N_orbit*T
N_out = 5
xy = np.zeros((N_out,N_testparticle,2))
times = np.linspace(0,t_max,N_out)
for i, time in enumerate(times):
	print i
	sim.integrate(time)
	for j, p in enumerate(sim.particles[8:]):
		print j
		xy[i][j] = [p.x, p.y]

#Now plot the positions of the test particles
#fig = plt.figure(figsize=(5,5))
#ax = plt.subplot(111)
plt.scatter(xy[:,:,0],xy[:,:,1],marker=".",linewidth='0')
#plt.xlim(-1e10,1e10)
#plt.ylim(-1e10,1e10)
plt.plot(0,0,'r.')
#plt.savefig('/home/james/Documents/URSS/Figures/TestParticles/Positions_'+str(planet)+'_'+str(N_orbit)+'orbits_ratio='+str(ratio)+'.png')
plt.show()

#Plot the relative changes in semimajor axis as a function of starting period
orbits = sim.calculate_orbits()[7:]
a_final = [o.a for o in orbits]
print len(a_final)
fig = plt.figure(figsize=(15,5))
ax = plt.subplot(111)
ax.set_yscale('log')
ax.set_xlabel(r"period ratio $r$")
ax.set_ylabel("relative semi-major axis change")
plt.plot(np.power(a_post,1.5),(np.fabs(a_final-a_post)+1.0e-16)/a_post,marker=".")
plt.show()
#Those showing a change of order unity or above have experienced a close encounter...
########################################################################
########################################################################
########################################################################
'''

########################################################################
##################### Look for close encounters ########################
########################################################################

#Calculate orbital period of planet for integration
T = orb_period(a)

#Integration set up
#sim.integrator = 'ias15' #Good for close encounters, although slower than hermes
sim.move_to_com()
sim.integrator = 'whfast'
sim.dt = 100

#Required to iterate through particle pairs
N_orbits = 1000
t_max = N_orbits*T
N_outputs = 1e6
times = np.linspace(0,t_max,N_outputs)
tp_indices = np.arange(sim.N_active,N_testparticle+sim.N_active,1)

d_star_init = {}
for i in tp_indices:
	d_star_init.update({i:[0]*1000})

merge_count = {}
def mergeParticles(sim, mappedParticles):
	ps = sim.particles
	for j, key in enumerate(sorted(pdict)):
		r_H = r_Hill[key]
		k = sim.N_active
		while k < sim.N:
			dp = ps[j] - ps[k]
			d2 = dp.x*dp.x+dp.y*dp.y+dp.z*dp.z
			d = np.sqrt(d2)
			
			if d < r_H:
				sim.remove(index=k)
				
				print 'Test particle %d has experienced a close encounter an merged with planet '%mappedParticles[k-sim.N_active][1]+key+'.'
				
				label = key
				merge_file(label,k)
				
				for l in range(k-sim.N_active, len(mappedParticles)-1):
					mappedParticles[l] = (mappedParticles[l][0], mappedParticles[l+1][1])
					
			if d > 1.5e11:
				sim.remove(index=k)
				
				print 'Test particle %d has experienced a long way from home.'%mappedParticles[k-sim.N_active][1]
				
				label = 'X'
				merge_file(label,k)
				
				for l in range(k-sim.N_active, len(mappedParticles)-1):
					mappedParticles[l] = (mappedParticles[l][0], mappedParticles[l+1][1])
			
			if key == 'A':
				d_star[mappedParticles[k-sim.N_active][1]][i % 1000] = d
			
			k += 1
			
mappedParticles = []
for i in tp_indices:
	mappedParticles.append((i,i))

d_star = d_star_init

keys = sorted(d_star.keys())
with open('/home/james/Documents/URSS/Output_files/'+system+'/d_star_'+planet+'_'+str(ratio)+'_'+str(N_testparticle)+'.csv', 'wb') as outfile:
	w = csv.writer(outfile, delimiter = '\t')
	w.writerow(keys)

with open('/home/james/Documents/URSS/Output_files/'+system+'/Merge_info_'+planet+'_'+str(ratio)+'_'+str(N_testparticle)+'.csv','a') as outfile:
	w = csv.writer(outfile, delimiter = '\t')
	w.writerow(['Index','Time','Outcome'])

def csv_writer(file_dict,keys):
    with open('/home/james/Documents/URSS/Output_files/'+system+'/d_star_'+planet+'_'+str(ratio)+'_'+str(N_testparticle)+'.csv', 'a') as outfile:
	    w = csv.writer(outfile, delimiter = '\t')
	    w.writerows(zip(*[d_star[key] for key in keys]))

def merge_file(label,k):
	with open('/home/james/Documents/URSS/Output_files/'+system+'/Merge_info_'+planet+'_'+str(ratio)+'_'+str(N_testparticle)+'.csv','a') as outfile:
		w = csv.writer(outfile,delimiter = '\t')
		w.writerow([mappedParticles[k-sim.N_active][1],t_orbits,label])

for i, time in enumerate(times):
    t = i * (t_max / N_outputs)
    t_orbits = t / T
    print t_orbits, sim.N
    
    if (i + 1) % 1000 == 0:
		csv_writer(d_star,keys)
		d_star = d_star_init
    
    sim.integrate(time)
    
    '''
    if i % 1000 == 0:
	    xy = np.zeros((N_testparticle+sim.N_active,2))   
	    for j, p in enumerate(sim.particles[:]):
		    xy[j] = [p.x, p.y]
	    
	    fig = plt.figure(figsize=(5,5))
	    ax = plt.subplot(111)
	    plt.scatter(xy[0,0],xy[0,1],marker="o",linewidth='0',color='r')
	    plt.scatter(xy[1:8,0],xy[1:8,1],marker="o",linewidth='0',color='k')
	    plt.scatter(xy[8:,0],xy[8:,1],marker=".",linewidth='0',color='b')
	    plt.xlim(-1e10,1e10)
	    plt.ylim(-1e10,1e10)
	    plt.savefig('/home/james/Documents/URSS/Figures/FirstStages/Positions_'+str(planet)+'_'+str(ratio)+'_'+str(N_testparticle)+'/'+str(t_orbits)+'.png')
	    #plt.show()
	    plt.close()
	'''
	
    if sim.N < sim.N_active + 1:
	    break
     
    mergeParticles(sim, mappedParticles)

########################################################################
########################################################################
########################################################################

'''
#Generalised functions for non-zero eccentricity - update later
def kep_vel():
	v_k = np.sqrt((G*(m_star + m))/a)
	return v_k

def ang_mom_ratio(phi):
    x = 1. + 2.*(np.sqrt(1 - np.square(e))/(1.+e*np.cos(f)))*ratio*np.sin(offset*(np.pi/180.))*np.sin(phi-f)
    y = ((1-np.square(e))/np.square((1+e*np.cos(f))))*np.square(ratio)*(np.square(np.cos(offset*(np.pi/180.))) + np.square(np.sin(offset*(np.pi/180.))*np.sin(phi-f)))
    return x + y

def post_incl(phi):
    H = ang_mom_ratio(phi)
    cos_I_post = (1. + (np.sqrt(1.-np.square(e))/(1+e*np.cos(f)))*ratio*np.sin(offset*(np.pi/180.))*np.sin(phi-f))*np.power(H,-0.5)
    I_post = np.arccos(cos_I_post)
    return I_post

def post_f(phi):
    num = np.abs(1. + ratio*np.sin(phi))*(1./np.tan(phi))
    den = 2. + ratio*np.sin(phi)
    tan_f_post = num / den
    f_post = np.arctan(tan_f_post)
    return f_post

def post_M(phi):
	f_post = post_f(phi)
	M_post = f_post - 2.*e*np.sin(f_post)
	return M_post

def post_a(phi):
	num = a
	den = 1. - np.square(ratio) - (2./np.sqrt(1-np.square(e)))*ratio*np.sin(offset*(np.pi/180.))*(np.sin(phi-f) + e*np.sin(phi))
	a_post = num / den
	return a_post

def post_e(phi):
	a_post = post_a(phi)
	H = ang_mom_ratio(phi)
	e_post = 1. - (1-np.square(e))*H*(a / a_post)
	return e_post
'''
