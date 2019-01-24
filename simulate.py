"""
Simulate the process of panspermia in an arbitrary exoplanetary system
"""

import rebound
import numpy as np
import argparse as ap

def argParse():
    """
    Argument parser settings
    
    Parameters
    ----------
    None
    
    Returns
    -------
    args : array-like
        Array of command line arguments
    """
    parser = ap.ArgumentParser()
    
    parser.add_argument('config',
                        help='path to config file',
                        type=str)
    
    parser.add_argument('time',
                        help='integration timescale [yrs]',
                        type=int)
    
    parser.add_argument('--diagnostics',
                        help='include sanity checks?',
                        action='store_true')
    
    return parser.parse_args()

if __name__ == "__main__":
    
    #args = argParse() 
    
    ####################################################################
    ############### Example from REBOUND documentation #################
    ####################################################################
    
    # setup simulation
    sim = rebound.Simulation()
    sim.integrator = "ias15" # not needed, ias15 is default integrator
    sim.add(m=1.)
    sim.add(m=1e-3,a=1.)
    sim.add(m=5e-3,a=1.25)
    sim.move_to_com()
    
    # can use exit_min_distance to quit upon close encounter
    # we'll instead want to track params as we go
    sim.exit_min_distance = 0.15
    Noutputs = 1000
    times = np.linspace(0,100.*2.*np.pi,Noutputs)
    distances = np.zeros(Noutputs)
    ps = sim.particles # ps is now an array of pointers 
                       # it will update as the simulation runs
    try:
        for i,time in enumerate(times):
            sim.integrate(time)
            dp = ps[1] - ps[2]   # calculates componentwise difference 
                                 # between particles
            distances[i] = np.sqrt(dp.x*dp.x+dp.y*dp.y+dp.z*dp.z)
        
    except rebound.Encounter as error:
        print(error)
    
    # plot distance over time
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(10,5))
    ax = plt.subplot(111)
    ax.set_xlabel("time [orbits]")
    ax.set_xlim([0,sim.t/(2.*np.pi)])
    ax.set_ylabel("distance")
    plt.plot(times/(2.*np.pi), distances);
    plt.plot([0.0,12],[0.2,0.2]); # plot our close encounter criteria
    
    plt.show()
    plt.close(fig)
