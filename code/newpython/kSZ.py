import numpy as np
import scipy as sp
import pylab as plt
from Simbox import *
from pspec import *
import timeit

def compute_tau(sim):
    """Computes the optical depth tau as a function of redshift for 
    a given simulation sim
    Input:  Sim object sim
    Output: 1D array tau(z)""" 

    # get proper distances and redshifts before the box begins
    z0,d0 = sim.redshift_to_space(0,sim.pms['zi'],sim.pms['shape'][2],
        proper=True)
    # get proper distances and redshifts for the box
    z1,d1 = sim.get_z_d(proper=True)
    # concatenate redshift and distances
    z = np.concatenate((z0,z1))
    d = np.concatenate((d0,d1))
    print 'have redshifts'

    # get the number density of baryons as a function of redshift
    nb = sim.pms['nb0']*(1+z)**3
    print 'have nb'
    
    # get the fraction of HII particles
    chiHII = np.ones_like(z)
    for kk in xrange(len(d)):
        start = timeit.default_timer()
        chiHII[kk+len(d0)] = np.average((1.-sim.box['nf'].slice(kk))* \
            (1.+sim.box['density'].slice(kk)))
        end = timeit.default_timer()
        print kk, end-start
    print 'have chiHII'

    # get the fraction of HeII particles
    chiHeIII = np.zeros_like(chiHII)
    chiHeIII[np.where(z<3)] = 1
    print 'have chiHeIII'

    # integrate to get tau(z)
    fn = sim.pms['sigT']*nb*(chiHII+
        0.25*chiHeIII*sim.pms['YBBNp'])*sim.pms['mToMpc']
    tau = sp.integrate.cumtrapz(fn,x=d,initial=0)
    print 'have tau'

    # only return tau(z) for the box
    return tau[len(d0):]


def compute_kSZ(sim):
    """Computes the kSZ temperature field.
    Input:  Sim object
    Output: 2D array of temperatures"""
    
    # get tau(z)
    tau = compute_tau(sim)
    
    # get z (redshift) and d
    zred,d = sim.get_z_d(proper=False)
    
    # do the integration slice by slice
    # import each slice as a numpy array during each iteration
    dTkSZ = np.zeros((sim.pms['shape'][0],sim.pms['shape'][1]))
    for kk in xrange(sim.pms['shape']):
        # Note: this version assumes line of sight is the z direction
        dTkSZ[:,:,kk] += np.exp(-tau[kk])*(1-sim.box['nf'].slice(kk))* \
            (1+sim.box['density'].slice(kk))* \
            sim.box['vz'].slice(kk)* \
            (1+zred[kk])*(d[kk+1]-d[kk])
    
    # Multiply by a bunch of constants
    dTkSZ = dTkSZ*box.sim.pms['sigT']*box.sim.pms['nb0']/box.sim.pms['c']*box.sim.pms['mToMpc']**2*box.sim.pms['Tcmb']
    return dTkSZ


if __name__=='__main__':
    data_dir = '/Users/mpresley/Research/KSZ_21cm_Constraints/data/mesinger_1/original_1_mesinger_1/'
    sim = Sim(data_dir) 
    tau = compute_tau(sim)
    z,d = sim.get_z_d()
    plt.plot(z,tau,c='orange')
    plt.xlabel('z')
    plt.ylabel('\\tau')
    plt.show()



    #dTkSZ = compute_kSZ(sim)




