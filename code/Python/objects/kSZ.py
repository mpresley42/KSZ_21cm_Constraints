import Sim21cm
from pspec import *
import numpy as np 
import scipy as sp
import pylab as plt 

# verified
def compute_tau(box):
    nf = box.get_data('nf')
    density = box.get_data('density')
    
    if box.ibox==0:
        z1,d1 = box.get_z_d(proper=True)
        z0,d0 = box.sim.get_z_d(0,z1[0],proper=True)
        z = np.concatenate((z0,z1))
        d = np.concatenate((d0,d1))
    else:
        z,d = box.get_z_d(proper=True)

    nb = box.sim.pms['nb0']*(1+z)**3
    chiHII = np.ones_like(z)
    if box.ibox==0: chiHII[len(d0):] = np.average((1.-nf)*(1+density),axis=(0,1))
    else:       chiHII = np.average((1.-nf)*(1+density),axis=(0,1))
    chiHeIII = np.zeros_like(chiHII)
    chiHeIII[np.where(z<3)] = 1

    fn = box.sim.pms['sigT']*nb*(chiHII+0.25*chiHeIII*box.sim.pms['YBBNp'])*box.sim.pms['mToMpc']
    tau = sp.integrate.cumtrapz(fn,x=d,initial=0)
    
    # NOTE: This returns the proper distance!
    if box.ibox==0: return z1,d1,tau[len(d0):]
    else:       return z,d,tau

# verified
def compute_ndotq(box,save=True):
    nf = box.get_data('nf')
    density = box.get_data('density')
    # nf, density
 # Find the density and ionization weighting for the velocity field q
    q0 = (1-nf)*(1+density)
    nf = None; density = None
    # q0
    z,d = box.get_z_d()
 # Get the spatial coordinates of all box cells
    xyzgd = box.get_xyz_gd()
    rgd = np.sqrt(xyzgd[0]*xyzgd[0]+xyzgd[1]*xyzgd[1]+xyzgd[2]*xyzgd[2])
    # q0, xyzgd (x3), rgd
 # Find the peculiar velocity 
    ndotq = np.zeros_like(q0)
    for ii,vii in enumerate(('vx','vy','vz')):
        vi = box.get_data(vii) #NOTE: These are co-moving velocities
        ndotq += q0*vi*xyzgd[ii]/rgd
    # q0, xyzgd (x3), rgd, ndotq, vi
    vi = None; q0=None; xyzgd=None; rgd=None
    if save: np.save('{0}ndotq_{1}'.format(box.sim.data_dir,box.ibox),ndotq)
    return ndotq

def compute_kSZ(box,ndotq=None,size=None,save=True):
    if ndotq==None: ndotq = compute_ndotq(box)
    print "Have ndotq!"
    # ndotq    
 # Define box shape parameters
    boxsh = box.ishape # shape of ibox
    boxMpc = np.array([box.sim.pms['xyMpc'],box.sim.pms['xyMpc'],box.sim.pms['zMpc']]) # size of ibox in Mpc
    lx = boxMpc[0]/2.; ly = boxMpc[1]/2.; lz = boxMpc[2]*box.itot # size of total box in Mpc
    dxyz = boxMpc/boxsh # length of each cell in Mpc
 # Get tau
    zred,dp,tau = compute_tau(box)
    zred,d = box.get_z_d(proper=False)
    # ndotq, tau
    # Get n,r coord
    nxcd,nycd,rcd = box.get_nxnyr_cd()
    if size==None: size = (len(nxcd),len(nycd))
    # Loop over angles
    dTkSZ = np.zeros([len(nxcd),len(nycd)])
    for ii in range(size[0]):
        print ii
        for jj in range(size[1]):
            # Loop over z direction
            for kk in range(len(rcd)-1):
                # Find corresponding xyz coord
                x = nxcd[ii]*rcd[kk]
                y = nycd[jj]*rcd[kk]
                z = np.sqrt(rcd[kk]*rcd[kk]-x*x-y*y)
                # Find closest xyz indices
                xi = int(np.floor((lx+x)/dxyz[0]))
                yj = int(np.floor((ly+y)/dxyz[1]))
                zk = int(np.floor((z-d[0])/dxyz[2]))
                if xi < boxsh[0] and yj < boxsh[1] and zk < boxsh[2]-1:
                    dTkSZ[ii,jj] += np.exp(-tau[zk])*ndotq[xi,yj,zk]*(1+zred[zk])*(d[zk+1]-d[zk])

    dTkSZ = dTkSZ*box.sim.pms['sigT']*box.sim.pms['nb0']/box.sim.pms['c']*box.sim.pms['mToMpc']**2*box.sim.pms['Tcmb']
    if save and size==(len(nxcd),len(nycd)): np.save('{0}dTkSZ_{1}'.format(box.sim.data_dir,box.ibox),dTkSZ)
    elif save: np.save('{0}dTkSZ_{1}_{2}_{3}'.format(box.sim.data_dir,size[0],size[1],box.ibox),dTkSZ)
    return dTkSZ 

def compute_kSZ_pspec(nxcd,nycd,dTkSZ,mask=None,pdw=512,n=100,pretty=True,save_dir='./'):
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
        plt.savefig('{0}dTkSZ.pdf'.format(save_dir))
        #plt.show()
 # Find the Fourier Transform
    lx,ly,dTkSZ_FT = fft_2d(nxcd,nycd,dTkSZ) # [K][Mpc]^2    
    if pretty:
        plt.clf()
        plt.imshow(np.real(dTkSZ_FT),origin='lower')
        cb=plt.colorbar();cb.set_label(r"$Real(\Delta \widetildeT_{kSZ})$",fontsize=18)
        plt.savefig('{0}dTkSZ_FT.pdf'.format(save_dir))
        #plt.show()
 # Find the Power Spectrum
    lbins,dTkSZ_P,area = pspec_2d(lx,ly,dTkSZ_FT,n=n) # [K]^2[Mpc]^4
    dTkSZ_P = dTkSZ_P/norm
    if pretty:
        plt.clf()
        plt.semilogy(lbins,dTkSZ_P)
        plt.ylabel(r"$P(\ell)$",fontsize=18)
        plt.xlabel(r"$\ell$",fontsize=18)
        plt.savefig('{0}dTkSZ_pspec.pdf'.format(save_dir))
        #plt.show()
    if pretty:
        plt.clf()
        #plt.semilogy(lbins,lbins*(lbins+1.)*dTkSZ_P/(2*np.pi)) # [Mpc]^-2[K]^2[Mpc]^4 = [K]^2[Mpc]^2
        plt.loglog(lbins,lbins*(lbins+1.)*dTkSZ_P/(2*np.pi))
        plt.ylabel(r"$\Delta_{kSZ}^2=\ell(\ell+1)C_\ell/2\pi [\mu K^2]$",fontsize=18)
        plt.xlabel(r"$\ell$",fontsize=18)
        plt.savefig('{0}dTkSZ_delta_kSZ.pdf'.format(save_dir))
        #plt.show()
    if True: np.save('{0}lbins_{1}'.format(save_dir,len(nxcd)),lbins)
    if True: np.save('{0}dTkSZ_P_{1}'.format(save_dir,len(nxcd)),dTkSZ_P)    
    return lbins,dTkSZ_P,area


