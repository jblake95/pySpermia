from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import csv
import numpy as np

system = str(raw_input('Planetary system: '))
planet = str(raw_input('Impacted planet: '))
N_testparticle = int(raw_input('Number of test particles: '))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for c, z in zip(['r', 'g', 'b'], [0.552, 0.276, 0.184]):
	
     times = []
     planets = []
     with open('C:\Users\James\Documents\Warwick\Year 4\URSS\Output_files\\'+system+'\\'+planet+'\Merge_info_'+planet+'_'+str(z)+'_'+str(N_testparticle)+'.csv','r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        for row in reader:
            times.append(row[1])
            planets.append(row[2])
        
        Targets = ['A','b','c','d','e','f','g','h','X']

        planet_counts = {}
    
        for i in Targets:
            planet_counts.update({i:planets.count(i)})
    
        counts = []
        for key in Targets:
            counts.append(planet_counts[key])
    
        remain = N_testparticle - len(times)

        Targets.append('-')
        counts.append(remain)
    
        Targets[8],Targets[9] = Targets[9],Targets[8]
        counts[8],counts[9] = counts[9],counts[8]
    
        y_pos = np.arange(len(Targets))
        
        if z == 0.184:
            z_pos = 1
        elif z == 0.276:
            z_pos = 2
        else:
            z_pos = 3
        
        cs = [c] * len(Targets)
    
        ax.bar(y_pos, counts, z_pos, zdir='y', color=cs, alpha=0.5)

z_labels = ['0.184','0.276','0.552']
z_pos = [1,2,3]
plt.xticks(y_pos, Targets)
#ax.set_ylim(0.5,2.5)
#ax.set_zlim(0,500)
ax.set_yticks(z_pos, z_labels)
ax.set_xlabel('Target')
ax.set_ylabel('Velocity ratio')
ax.set_zlabel('% test particles')
plt.show()

		
