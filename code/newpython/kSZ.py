import numpy as np
import scipy as sp
import pylab as plt
import Simbox
from pspec import *

def compute_tau(sim):
    """Computes the optical depth tau as a function of redshift for 
    a given simulation sim
    Input:  Sim object sim
    Output: 1D array tau(z)""" 

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
        dTkSZ[:,:,kk] += np.exp(-tau[kk])*(1-sim.box['nf'].slice(kk))* \
            (1+sim.box['density'].slice(kk))* \
            # Note: this version assumes line of sight is the z direction
            sim.box['vz'].slice(kk)* \
            (1+zred[kk])*(d[kk+1]-d[kk])
    
    # Multiply by a bunch of constants
    dTkSZ = dTkSZ*box.sim.pms['sigT']*box.sim.pms['nb0']/box.sim.pms['c']*box.sim.pms['mToMpc']**2*box.sim.pms['Tcmb']
    return dTkSZ