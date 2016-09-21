import sys
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
    for ii in range(sim.itot):
        box = Box21cm(sim,ii)
        ndotq = np.load('{0}ndotq_{1}.npy'.format(sim.data_dir,ii))
        kSZ.compute_kSZ(box)
        print "computed dTkSZ for box {0}".format(ii)

def compute_full_kSZ_pspec(sim):
    dTkSZ = np.zeros((sim.pms['ishape'][0],sim.pms['ishape'][0]))
    for ii in range(sim.itot):
        dTkSZ += np.load('{0}dTkSZ_256_256_{1}.npy'.format(sim.data_dir,ii))
    box = Box21cm(sim,0)
    nxcd,nycd,rcd = box.get_nxnyr_cd()
    dTkSZ_P = kSZ.compute_kSZ_pspec(nxcd[:256],nycd[:256],dTkSZ[:256,:256],save_dir=sim.data_dir+'figures/')

if __name__=='__main__':
    #ibox = int(sys.argv[1])
    sim = Sim21cm(data_dir,run_dir,itot=8)
    compute_full_kSZ_pspec(sim)
    # box = Box21cm(sim,ibox)
    # ndotq = np.load('{0}ndotq_{1}.npy'.format(sim.data_dir,ibox))
    # kSZ.compute_kSZ(box,ndotq,size=(256,256))
    
    #dTkSZ = np.load('{0}dTkSZ_256_256_{0}.npy'.format(ibox)))


