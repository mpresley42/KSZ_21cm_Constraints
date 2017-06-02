import numpy as np
import scipy as sp
import pylab as plt
from Simbox import *
from pspec import *
import timeit
import sys

def scrape_data(sim):
    """Collects z, nf, and Tb data from the filenames of 21cmFAST
    data boxes and writes them to a .dat file."""
    pattern = 'delta_T_v3_*Mpc'
    filenames = match_files(sim.data_dir,pattern)
    wf = open('{0}z_nf_Tb.dat'.format(location),'w')
    for f in filenames:
        nums = nums_from_string(f)
        z = nums[1]
        nf = nums[2]
        Tb = nums[9]
        wf.write('{0} {1} {2}\n'.format(z,nf,Tb))
    wf.close()

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

    # get the number density of baryons as a function of redshift
    nb = sim.pms['nb0']*(1+z)**3
    #print 'have nb'
    
    # get the fraction of HII particles
    chiHII = np.ones_like(z)
    for kk in xrange(len(d1)):
        chiHII[kk+len(d0)] = np.average((1.-sim.box['nf'].slice(kk))* \
            (1.+sim.box['density'].slice(kk)))

    # get the fraction of HeII particles
    chiHeIII = np.zeros_like(chiHII)
    chiHeIII[np.where(z<3)] = 1

    # integrate to get tau(z)
    fn = sim.pms['sigT']*nb*(chiHII+
        0.25*chiHeIII*sim.pms['YBBNp'])*sim.pms['mToMpc']
    tau = sp.integrate.cumtrapz(fn,x=d,initial=0)
    # select only the tau(z) inside the box
    tau = tau[len(d0):]

    # save and return tau
    np.save('{0}tau'.format(sim.data_dir),tau)
    return tau

def compute_kSZ(sim,tau=None):
    """Computes the kSZ temperature field.
    Input:  Sim object
    Output: 2D array of temperatures"""
    
    # get tau(z)
    if tau is None: tau = compute_tau(sim)
    
    # get z (redshift) and d
    zred,d = sim.get_z_d(proper=False)
    
    # do the integration slice by slice
    # import each slice as a numpy array during each iteration
    start = timeit.default_timer()
    dTkSZ = np.zeros((sim.pms['shape'][0],sim.pms['shape'][1]))
    for kk in xrange(sim.pms['shape'][2]-1):
        # Note: this version assumes line of sight is the z direction
        dTkSZ += np.exp(-tau[kk])*(1-sim.box['nf'].slice(kk))* \
            (1+sim.box['density'].slice(kk))* \
            sim.box['vz'].slice(kk)* \
            (1+zred[kk])*(d[kk+1]-d[kk])
    end = timeit.default_timer()
    print 'dTkSZ integration: ', end-start
    
    # Multiply by a bunch of constants
    dTkSZ = dTkSZ*sim.pms['sigT']*sim.pms['nb0']/sim.pms['c']*sim.pms['mToMpc']**2*sim.pms['Tcmb']

    # save and return dTkSZ field
    np.save('{0}dTkSZ'.format(sim.data_dir),dTkSZ)
    return dTkSZ


def compute_kSZ_pspec(sim,dTkSZ=None,mask=None,pdw=512,n=100,pretty=True):
    """Computes the kSZ power spectrum from the temperature field.
    Input:  Sim object 
    Output: 1D array of lbins, 1D array of power per bin"""
 
    # get dTkSZ
    if dTkSZ is None: dTkSZ = compute_kSZ(sim)

    # get nx,ny coordinates
    nxcd,nycd,_ = sim.get_nxnyr_cd()

    # Apply mask if necessary 
    if mask==None: mask = np.ones_like(dTkSZ)
    else: dTkSZ = dTkSZ*mask
    
    # Compute overall normalization
    norm = np.trapz(np.trapz(mask*mask,nycd,axis=1),nxcd,axis=0)
    
    # Apply padding if necessary
    nxcd,nycd,dTkSZ = pad2_array(nxcd,nycd,dTkSZ,pdw)
    
    # Plot the new dTkSZ 
    if pretty:
        plt.clf()
        plt.imshow(dTkSZ,origin='lower')
        cb=plt.colorbar();cb.set_label(r"$\Delta T_{kSZ}$",fontsize=18)
        plt.savefig('{0}dTkSZ.pdf'.format(sim.data_dir))
    
    # Find the Fourier transform
    lx,ly,dTkSZ_FT = fft_2d(nxcd,nycd,dTkSZ) # [K][Mpc]^2    
    if pretty:
        plt.clf()
        plt.imshow(np.real(dTkSZ_FT),origin='lower')
        cb=plt.colorbar();cb.set_label(r"$Real(\Delta \widetildeT_{kSZ})$",fontsize=18)
        plt.savefig('{0}dTkSZ_FT.pdf'.format(sim.data_dir))

    # Find the power spectrum
    lbins,dTkSZ_P,area = pspec_2d(lx,ly,dTkSZ_FT,n=n) # [K]^2[Mpc]^4
    dTkSZ_P = dTkSZ_P/norm

    # Plot the power spectrum
    if pretty:
        plt.clf()
        plt.semilogy(lbins,dTkSZ_P)
        plt.ylabel(r"$P(\ell)$",fontsize=18)
        plt.xlabel(r"$\ell$",fontsize=18)
        plt.savefig('{0}dTkSZ_pspec.pdf'.format(sim.data_dir))

    # Find Delta-squared (the renormalized power spectrum)
    delta_sq = lbins*(lbins+1.)*dTkSZ_P/(2*np.pi)
    
    # Plot Delta-squared
    if pretty:
        plt.clf()
        #plt.semilogy(lbins,delta_sq) # [Mpc]^-2[K]^2[Mpc]^4 = [K]^2[Mpc]^2
        plt.loglog(lbins,delta_sq)
        plt.ylabel(r"$\Delta_{kSZ}^2=\ell(\ell+1)C_\ell/2\pi [\mu K^2]$",fontsize=18)
        plt.xlabel(r"$\ell$",fontsize=18)
        plt.savefig('{0}dTkSZ_delta_kSZ.pdf'.format(sim.data_dir))

    # Save the lbins, power spectrum, and Delta-squared
    np.save('{0}lbins'.format(sim.data_dir),lbins)
    np.save('{0}dTkSZ_P'.format(sim.data_dir),dTkSZ_P)  
    np.save('{0}Delta_sq'.format(sim.data_dir),delta_sq)  
    return lbins,dTkSZ_P

def timing_test():
    """Test how long it takes to compute different steps of the 
    kSZ calculation."""
    data_dir = '/Users/mpresley/Research/KSZ_21cm_Constraints/data/small/RandSeed_111_Sigma8_0.81577_h_0.68123_Omm_0.30404_Omb_0.04805_ns_0.96670_Rmfp_35.00_Tvir_60000.0_Squiggly_40.00_lnAs_3.06400/'

    tot_start = timeit.default_timer()
    # initialize simbox
    start = timeit.default_timer()
    sim = Sim(data_dir) 
    end = timeit.default_timer()
    print 'Simbox initialization: ', end-start
    
    # get tau
    start = timeit.default_timer()
    tau = compute_tau(sim)    
    np.save('tau.npy',tau)
    end = timeit.default_timer()
    print 'tau computation: ', end-start
    #tau = np.load('tau.npy')

    # ksz calc
    start = timeit.default_timer()
    dTkSZ = compute_kSZ(sim,tau)
    end = timeit.default_timer()
    print 'dTkSZ computation: ', end-start

    # pspec calc
    start = timeit.default_timer()
    dTkSZ_P = compute_kSZ_pspec(sim,dTkSZ)
    end = timeit.default_timer()
    print 'Pspec computation: ', end-start

    tot_end = timeit.default_timer()
    print 'total time: ',tot_end - tot_start

if __name__=='__main__':
    #data_dir = '/Users/mpresley/Research/KSZ_21cm_Constraints/data/mesinger_1/original_1_mesinger_1/'
    #data_dir = '/Users/mpresley/Research/KSZ_21cm_Constraints/data/small/RandSeed_111_Sigma8_0.81577_h_0.68123_Omm_0.30404_Omb_0.04805_ns_0.96670_Rmfp_35.00_Tvir_60000.0_Squiggly_40.00_lnAs_3.06400/'

    if len(sys.argv)==1:
        data_dir = './'
    else:
        data_dir = sys.argv[1]    

    sim = Sim(data_dir)
    scrape_data(sim)
    compute_kSZ_pspec(sim)
    



