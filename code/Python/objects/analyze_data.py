import numpy as np 
import pylab as plt
from Sim21cm import *
import kSZ 

data_dir = '/Users/mpresley/Research/KSZ_21cm_Constraints/data/mesinger/'
run_dir = 'mesinger_1'

def plot_tau(sim):
    tau0 = 0.
    for ii in range(8):
        box = Box21cm(sim,ii)
        z,dp,tau = kSZ.compute_tau(box)
        plt.plot(z,tau+tau0)
        tau0 = tau[-1]+tau0
    plt.xlabel(r'$z$',fontsize=18)
    plt.ylabel(r'$\tau$',fontsize=18)
    plt.show()

def compute_all_ndotq(sim):
    for ii in range(2,8):
        box = Box21cm(sim,ii)
        kSZ.compute_ndotq(box)
        print "computed ndotq for box {0}".format(ii)

def compute_all_kSZ(sim):
    for ii in range(8):
        box = Box21cm(sim,ii)
        ndotq = np.load('{0}ndotq_{1}.npy'.format(sim.data_dir,ii))
        kSZ.compute_kSZ(box)
        print "computed dTkSZ for box {0}".format(ii)


if __name__=='__main__':
    sim = Sim21cm(data_dir,run_dir,itot=8)
    box = Box21cm(sim,0)
    ndotq = np.load('{0}ndotq_0.npy'.format(sim.data_dir))
    kSZ.compute_kSZ(box,ndotq,size=(256,256))

