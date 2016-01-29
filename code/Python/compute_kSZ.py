import cfg
import numpy as np 
import scipy as sp
import scipy.integrate
from load_data import *
from redshift import *
from scipy.interpolate import griddata
from warnings import warn 

def compute_tau(density,nf):
    z,d = get_z_d(cfg.pms['zi'],cfg.pms['zf'])
    nb = cfg.pms['nb0']*(1+z)**3
    chiHII = np.average(nf*(1+density),axis=(0,1))
    chiHeIII = np.zeros_like(chiHII)
    chiHeIII[np.where(z<3)] = 1
    fn = cfg.pms['sigT']*nb*(chiHII+0.25*chiHeIII*cfg.pms['YBBNp'])*mToMpc
    tau = sp.integrate.cumtrapz(fn,dx=(d[-1]-d[0])/len(d),initial=0)
    return z,d,tau

def compute_ndotq():
    nf = get_int_box_data('nf')
    density = get_int_box_data('density')
    # nf, density
 # Find the density and ionization weighting for the velocity field q
    q0 = (1+nf)*(1+density)
    nf = None; density = None
    # q0
    z,d = get_z_d(cfg.pms['zi'],cfg.pms['zf'])
 # Get the spatial coordinates of all box cells
    xyzgd = get_xyz_gd()
    rgd = np.sqrt(xyzgd[0]*xyzgd[0]+xyzgd[1]*xyzgd[1]+xyzgd[2]*xyzgd[2])
    # q0, xyzgd (x3), rgd
 # Find the peculiar velocity 
    ndotq = np.zeros_like(q0)
    for ii,vii in enumerate(('vx','vy','vz')):
        vi = get_int_box_data(vii)
        ndotq += q0*vi*xyzgd[ii]/rgd
    # q0, xyzgd (x3), rgd, ndotq, vi
    vi = None; q0=None; xyzgd=None; rgd=None
    if True: np.save('ndotq',ndotq)
    return ndotq

# old kSZ code
def compute_kSZ():
    nf = get_int_box_data('nf')
    density = get_int_box_data('density')
    # nf, density
 # Find the density and ionization weighting for the velocity field q
    q0 = (1+nf)*(1+density)
    nf = None; density = None
    # q0
 # Get the spatial coordinates of all box cells
    z,d = get_z_d(cfg.pms['zi'],cfg.pms['zf'])
    box = q0.shape
    xcd = np.linspace(-box[0]/2,box[0]/2,num=box[0])
    ycd = np.linspace(-box[1]/2,box[1]/2,num=box[1])
    zcd = d[0] + np.linspace(-box[2]/2,box[2]/2,num=box[2])
    xyzgd = np.meshgrid(xcd,ycd,zcd)
    xcd,ycd,zcd = None,None,None
    rgd = np.sqrt(xyzgd[0]*xyzgd[0]+xyzgd[1]*xyzgd[1]+xyzgd[2]*xyzgd[2])
    # q0, xyzgd (x3), rgd
 # Find the peculiar velocity 
    ndotq = np.zeros_like(q0)
    for ii,vii in enumerate(('vx','vy','vz')):
        vi = get_int_box_data(vii)
        ndotq += q0*vi*xyzgd[ii]/rgd
    # q0, xyzgd (x3), rgd, ndotq, vi
    vi = None; q0=None; rgd=None; xyzgd=None;
    # ndotq
 # Get tau
    nf = get_int_box_data('nf')
    density = get_int_box_data('density')
    z,d,tau = compute_tau(density,nf)
    # ndotq, nf, density, tau
    nf=None; density=None
    # ndotq, tau
 # Compute the kSZ signal
    fn = np.exp(-tau)*ndotq/cfg.pms['c']*mToMpc*cfg.pms['Tcmb']*(1+z) # (Mpc/s)(m/s)^-1(m/Mpc)(K) = K
    dTkSZ = sp.integrate.trapz(fn,x=tau,axis=2) # K
    fn = None
    return dTkSZ 

def get_xyz_gd():
    box = cfg.pms['shape']
    xcd = np.linspace(-box[0]/2,box[0]/2,num=box[0])
    ycd = np.linspace(-box[1]/2,box[1]/2,num=box[1])
    z,d = get_z_d(cfg.pms['zi'],cfg.pms['zf'])
    zcd = d[0] + np.linspace(-box[2]/2,box[2]/2,num=box[2])
    xyzgd = np.meshgrid(xcd,ycd,zcd)
    xcd,ycd,zcd = None,None,None
    return xyzgd

def get_nxnyr_gd():
    xgd,ygd,zgd = get_xyz_gd()
    rgd = np.sqrt(xgd*xgd+ygd*ygd+zgd*zgd)
    nygd = ygd/rgd
    nycd = np.linspace(np.amin(nygd),np.amax(nygd),cfg.pms['shape'][1])
    nygd = None; ygd=None
    nxgd = xgd/rgd
    nxcd = np.linspace(np.amin(nxgd),np.amax(nxgd),cfg.pms['shape'][0])
    nxgd = None; xgd=None
    rcd = np.linspace(np.amin(rgd),np.amax(rgd),cfg.pms['shape'][2])
    rgd = None
    ngd = np.meshgrid(nxcd,nycd,rcd)
    return ngd

# Function to find indices from a position
def invert_linspace(a,b,n,x,c=0):
    """a,b,n are the inputs to linspace. x is the point you want to invert.
       c is a contant offset. the index i is returned.
       x = a + (b-a)i/(n-1) + c ==> i = (n-1)(x-a-c)/(b-a)"""
    i = int(round((n-1)*float(x-a-c)/(b-a)))
    if i < 0: 
        warn("Index too small (i = {0}).".format(i))
        return 0
    elif i >= n: 
        warn("Index too big (i = {0}).".format(i))
        return n-1
    else: 
        return i
def get_xyz_id(xp,box,cs=(0,0,0)):
    xid = invert_linspace(-box[0]/2,box[0]/2,box[0],xp[0],cs[0])
    yid = invert_linspace(-box[1]/2,box[1]/2,box[1],xp[1],cs[1])
    zid = invert_linspace(-box[2]/2,box[2]/2,box[2],xp[2],cs[2])
    return (xid,yid,zid)

def regrid_old(ndotq):
 # create new nx,ny,r coordinate grids
    ngd = get_nxnyr_gd()
 # find the xyz values that correspond to the nx,ny,r coordinate grid
    n_xgd = ngd[0]*ngd[2]
    n_ygd = ngd[1]*ngd[2]
    n_zgd = np.sqrt(ngd[2]*ngd[2] - ngd[0]*ngd[0] - ngd[1]*ngd[1])
    ngd=None
 # create the new grid of ndotq
    ndotq_ncd = np.zeros(cfg.pms['shape'])
    print "Starting loop!"
    for ind in np.ndindex(ndotq.shape):
        xyz_id = get_xyz_id((n_xgd[ind],n_ygd[ind],n_zgd[ind]),cfg.pms['shape'])
        if ind[0]%10==0 and ind[1]==0 and ind[2]==0: 
            print ind, xyz_id
        ndotq_ncd[ind] = ndotq[xyz_id]
    # loop runs for 17s on a grid of 10x400x400
    # should run for 11.3 min on a 400x400x400 grid
    if True: np.save('ndotq_ncd',ndotq_ncd)
    return ndotq_ncd

def regrid(ndotq):
 # create new nx,ny,r coordinate grids
    ngd = get_nxnyr_gd()
 # find the xyz values that correspond to the nx,ny,r coordinate grid
    n_xgd = ngd[0]*ngd[2]
    n_ygd = ngd[1]*ngd[2]
    #n_zgd = np.sqrt(ngd[2]*ngd[2] - ngd[0]*ngd[0] - ngd[1]*ngd[1])
    n_zgd = np.sqrt(ngd[2]*ngd[2] - n_xgd*n_xgd - n_ygd*n_ygd)
    ngd=None
 # find the indices of the closest point on the regular xyz grid
    z,d = get_z_d(cfg.pms['zi'],cfg.pms['zf'])
    box = np.array(cfg.pms['shape'])
    xyz_0 = -box/2.+np.array([0,0,d[0]])
    xyz_n = box/2.+np.array([0,0,d[0]])
    xyz_delta = (xyz_n - xyz_0)/box
    n_igd = ((n_xgd - xyz_0[0])/xyz_delta[0]).round().clip(0,box[0]-1).astype(int)
    n_jgd = ((n_ygd - xyz_0[1])/xyz_delta[1]).round().clip(0,box[1]-1).astype(int)
    n_kgd = ((n_zgd - xyz_0[2])/xyz_delta[2]).round().clip(0,box[2]-1).astype(int)
    print np.amax(n_igd),np.amax(n_jgd),np.amax(n_kgd)
    print np.amin(n_igd),np.amin(n_jgd),np.amin(n_kgd)
    # plt.imshow(n_igd[:,:,0],origin='lower'); plt.show()
    # plt.imshow(n_jgd[:,:,0],origin='lower'); plt.show()
    # plt.imshow(n_kgd[0,:,:],origin='lower'); plt.show()
 # create the new grid of ndotq
    ndotq_ncd = ndotq[n_igd,n_jgd,n_kgd]
    if True: np.save('ndotq_ncd',ndotq_ncd)
    return ndotq_ncd

def compute_kSZ2(ndotq=None):
    if ndotq==None: ndotq = compute_ndotq()
    print "Have ndotq!"
    # ndotq    
 # Get tau
    nf = get_int_box_data('nf')
    density = get_int_box_data('density')
    z,d,tau = compute_tau(density,nf)
    # ndotq, nf, density, tau
    nf=None; density=None
    # ndotq, tau
    print "Beginning regrid!"
    ndotq_ncd = regrid(ndotq)
    ndotq=None
    print "Regrid complete!"
    # ndotq_ncd, tau
 # Compute the kSZ signal
    fn = np.exp(-tau)*ndotq_ncd/cfg.pms['c']*mToMpc*cfg.pms['Tcmb']*(1+z) # (Mpc/s)(m/s)^-1(m/Mpc)(K) = K
    dTkSZ = sp.integrate.trapz(fn,x=tau,axis=2) # K
    # ^ an array w/ axes theta and phi
    fn = None
    return dTkSZ 

# 2d Fourier Transform (note: this uses the numpy convention, not the cosmology convention)
# Units: [k]  = 1/[x]
#        [ft] = [fg][x]^2 = [fg][k]^-2
def fft_2d(x,y,fg):
    lx = max(x)-min(x); ly = max(y)-min(y)
    nx = float(len(x)); ny = float(len(y))
    ftd = np.fft.fftshift(np.fft.fftn(np.fft.fftshift(fg)))*(lx/nx)*(ly/ny)
    kx = 2.*np.pi*np.fft.fftshift(np.fft.fftfreq(int(nx),d=lx/nx))
    ky = 2.*np.pi*np.fft.fftshift(np.fft.fftfreq(int(ny),d=ly/ny))
    return kx,ky,ftd

# Compute Power Spectrum from a Fourier Transform
# Units: [kbins] = [k]
#        [pspec] = [ft]^2 = [fg]^2[x]^4
def pspec_2d(kx,ky,ft,n=100):
    ft = np.abs(ft)**2
    kxg,kyg = np.meshgrid(kx,ky)
    kg = np.sqrt(kxg*kxg+kyg*kyg)
    kmin = np.min(kg);kmax = np.max(kg)
    kbins = np.arange(kmin,kmax,(kmax-kmin)/n)
    pspec = np.zeros_like(kbins)
    for ii in range(len(kbins)-1):
        pspec[ii] = np.sum(np.where(np.logical_and(kg>=kbins[ii], kg<kbins[ii+1]), ft, 0))
        area = np.sum(np.where(np.logical_and(kg>=kbins[ii], kg<kbins[ii+1]), 1, 0))
        pspec[ii] = pspec[ii]/area
    kbins=kbins[:-1]; pspec=pspec[:-1]
    return kbins,pspec

def compute_kSZ_pspec(dTkSZ):
    #x = np.arange(cfg.pms['xyMpc']) - cfg.pms['xyMpc']/2
 # Get the spatial xyz coords
    z,d = get_z_d(cfg.pms['zi'],cfg.pms['zf'])
    box = cfg.pms['shape']
    x = np.linspace(-box[0]/2,box[0]/2,num=box[0])
    y = np.linspace(-box[1]/2,box[1]/2,num=box[1])
    z = d[0] + np.linspace(-box[2]/2,box[2]/2,num=box[2])
 # Get the nx, ny coords
    R = np.sqrt(x*x + y*y + z*z)
    nx = x/R
    ny = y/R
 # Find the Fourier Transform
    kx,ky,dTkSZ_FT = fft_2d(x,x,dTkSZ) # [K][Mpc]^2
    dTkSZ_FT = np.abs(dTkSZ_FT)    
    if False:
        print dTkSZ_FT.min(), dTkSZ_FT.max(), dTkSZ_FT.mean()
        plt.imshow(np.log(dTkSZ_FT),cmap='Blues',origin='lower')
        cb=plt.colorbar();cb.set_label(r"$Log(\Delta \widetildeT_{kSZ})$",fontsize=18)
        plt.show()
 # Find the Power Spectrum
    kbins,dTkSZ_P = pspec_2d(kx,ky,dTkSZ_FT) # [K]^2[Mpc]^4
    if False:
        plt.scatter(kbins,np.log(dTkSZ_P))
        plt.ylabel(r"$\log[P(k)]$",fontsize=18)
        plt.xlabel(r"$k$",fontsize=18)
        plt.show()
    if True:
        plt.scatter(kbins,kbins*kbins*dTkSZ_P/(2.*np.pi)) # [Mpc]^-2[K]^2[Mpc]^4 = [K]^2[Mpc]^2
        plt.ylabel(r"$\Delta_{kSZ}^2=\ell(\ell+1)C_\ell/2\pi [\mu K^2]$",fontsize=18)
        plt.xlabel(r"$k$",fontsize=18)
        plt.show()

if __name__=='__main__':
    # test for compute_tau()
    nf = get_int_box_data('nf')
    # density = get_int_box_data('density')
    # z,d,tau = compute_tau(density,nf)
    # plt.plot(z,tau); plt.show()

    # test for compute_kSZ()
    ndotq = np.load('ndotq.npy')
    dTkSZ = compute_kSZ2(ndotq)
    #dTkSZ = compute_kSZ2()
    plt.imshow(dTkSZ,origin='lower')
    plt.xlabel(r"$x\ (\mathrm{Mpc})$",fontsize=18)
    plt.ylabel(r"$y\ (\mathrm{Mpc})$",fontsize=18)
    cb=plt.colorbar()
    cb.set_label(r"$\Delta T_{kSZ}$",fontsize=18)
    plt.show()
    np.save('kSZ_array2',dTkSZ)

    # dTkSZ = np.load('kSZ_array.npy')
    # compute_kSZ_pspec(dTkSZ)

