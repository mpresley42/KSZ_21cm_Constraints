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

def compute_full_kSZ_pspec(sim,size=(256,256),linear=True):
    dTkSZ = np.zeros((sim.pms['ishape'][0],sim.pms['ishape'][0]))
    for ii in range(sim.itot):
        if linear: dTkSZ += np.load('{0}dTkSZ_linear_{1}_{2}_{3}.npy'.format(sim.data_dir,size[0],size[1],ii))
        else:      dTkSZ += np.load('{0}dTkSZ_{1}_{2}_{3}.npy'.format(sim.data_dir,size[0],size[1],ii))
    box = Box21cm(sim,0)
    # nxcd,nycd,rcd = box.get_nxnyr_cd()
    # lbins,dTkSZ_P = kSZ.compute_kSZ_pspec(nxcd[:size[0]],nycd[:size[1]],dTkSZ[:size[0],:size[1]],save_dir=sim.data_dir+'figures/linear/',pretty=False)
    xcd,ycd,zcd = box.get_xyz_cd()
    lbins,dTkSZ_P = kSZ.compute_kSZ_pspec(xcd[:size[0]],ycd[:size[1]],dTkSZ[:size[0],:size[1]],save_dir=sim.data_dir+'figures/linear/',pretty=False)
    return lbins,dTkSZ_P

def compare_linear_ray(sim):
    lbins,dTkSZ_P = compute_full_kSZ_pspec(sim,size=(256,256),linear=True)
    plt.loglog(lbins,lbins*(lbins+1.)*dTkSZ_P/(2*np.pi),label='z integration') 
    lbins,dTkSZ_P = compute_full_kSZ_pspec(sim,size=(256,256),linear=False)
    plt.loglog(lbins,lbins*(lbins+1.)*dTkSZ_P/(2*np.pi),label='ray integration')

    plt.ylabel(r"$\Delta_{kSZ}^2=\ell(\ell+1)C_\ell/2\pi [\mu K^2]$",fontsize=18)
    plt.xlabel(r"$\ell$",fontsize=18)
    plt.legend()
    plt.savefig('{0}dTkSZ_delta_kSZ_comparison.pdf'.format(sim.data_dir+'figures/linear/'))

def test_ndotq_linear(sim):
    box = Box21cm(sim,0)
    nf = box.get_data('nf')
    density = box.get_data('density')
    vz = box.get_data('vz')
    ndotq = (1-nf)*(1+density)*vz
    nf = None; density = None; vz = None
    np.save('{0}ndotq_linear_{1}'.format(box.sim.data_dir,box.ibox),ndotq)


if __name__=='__main__':
    # ibox = int(sys.argv[1])
    # sim = Sim21cm(data_dir,run_dir,itot=8)
    # box = Box21cm(sim,ibox)
    # ndotq = np.load('{0}ndotq_{1}.npy'.format(sim.data_dir,ibox))
    # kSZ.compute_kSZ(box,ndotq,size=(256,256),save=True)

    sim = Sim21cm(data_dir,run_dir,itot=8)
    compare_linear_ray(sim)
    
    
    # ksz=0.;ksz_lin=0.
    # for ibox in range(0,1):#sim.itot):
    #     print ibox
    #     box = Box21cm(sim,ibox)
    #     ndotq = np.load('{0}ndotq_{1}.npy'.format(sim.data_dir,ibox))
    #     ksz += kSZ.compute_kSZ(box,ndotq=ndotq,size=None,save=False)[512,512]
    #     ksz_lin += np.load('{0}dTkSZ_linear_1024_1024_{1}.npy'.format(sim.data_dir,ibox))[512,512]
    # print 'ksz_ray = {0}'.format(ksz)
    # print 'ksz_lin = {0}'.format(ksz_lin)

    #compute_full_kSZ_pspec(sim,size=size)

