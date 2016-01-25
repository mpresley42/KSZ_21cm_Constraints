import cfg
import numpy as np 
import scipy as sp
import scipy.integrate
from load_data import *
from redshift import *
from scipy.interpolate import griddata

def compute_tau(density,nf):
    z,d = get_z_d(cfg.pms['zi'],cfg.pms['zf'])
    nb = cfg.pms['nb0']*(1+z)**3
    chiHII = np.average(nf*(1+density),axis=(0,1))
    chiHeIII = np.zeros_like(chiHII)
    chiHeIII[np.where(z<3)] = 1
    fn = cfg.pms['sigT']*nb*(chiHII+0.25*chiHeIII*cfg.pms['YBBNp'])*mToMpc
    tau = sp.integrate.cumtrapz(fn,dx=(d[-1]-d[0])/len(d),initial=0)
    return z,d,tau

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

def compute_kSZ2():
    nf = get_int_box_data('nf')
    density = get_int_box_data('density')
    # nf, density
 # Find the density and ionization weighting for the velocity field q
    q0 = (1+nf)*(1+density)
    nf = None; density = None
    # q0
    z,d = get_z_d(cfg.pms['zi'],cfg.pms['zf'])
 # Get the spatial coordinates of all box cells
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
    vi = None; q0=None; 
    # xyzgd (x3), rgd, ndotq    
 # Get tau
    nf = get_int_box_data('nf')
    density = get_int_box_data('density')
    z,d,tau = compute_tau(density,nf)
    # xyzgd (x3), ndotq, nf, density, tau
    nf=None; density=None
    # xyzgd (x3), ndotq, tau
 # Convert to spherical coords
    thgd = np.arccos(xyzgd[2]/rgd)
    phgd = np.arctan(xyzgd[1]/xyzgd[0])
    xyzgd=None;
    # rgd, thgd, phgd, ndotq, tau
    thcd = np.linspace(np.amin(thgd),np.amax(thgd),cfg.pms['shape'][0])
    phcd = np.linspace(np.amin(phgd),np.amax(phgd),cfg.pms['shape'][1])
    rcd = np.linspace(np.amin(rgd),np.amax(rgd),cfg.pms['shape'][2])
    ndotq_s = griddata((thgd, phgd, rgd), ndotq, (thcd, phcd, rcd), method='linear')
    ndotq=None;
    # rgd, thgd, phgd, ndotq_s, tau
 # Compute the kSZ signal
    fn = np.exp(-tau)*ndotq_s/cfg.pms['c']*mToMpc*cfg.pms['Tcmb']*(1+z) # (Mpc/s)(m/s)^-1(m/Mpc)(K) = K
    dTkSZ = sp.integrate.trapz(fn,x=tau,axis=2) # K
    # ^ an array w/ axes theta and phi
    fn = None
    return dTkSZ 

def compute_kSZ3():
    nf = get_int_box_data('nf')
    density = get_int_box_data('density')
    # nf, density
 # Find the density and ionization weighting for the velocity field q
    q0 = (1+nf)*(1+density)
    nf = None; density = None
    # q0
    z,d = get_z_d(cfg.pms['zi'],cfg.pms['zf'])
 # Get the spatial coordinates of all box cells
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
    vi = None; q0=None; 
    # xyzgd (x3), rgd, ndotq    
 # Get tau
    nf = get_int_box_data('nf')
    density = get_int_box_data('density')
    z,d,tau = compute_tau(density,nf)
    # xyzgd (x3), ndotq, nf, density, tau
    nf=None; density=None
    # xyzgd (x3), ndotq, tau
 # Convert to nx,ny,r coords
    nxgd = xyzgd[0]/rgd
    nygd = xyzgd[1]/xyzgd[0]
    xyzgd=None;
    # rgd, nxgd, nygd, ndotq, tau
    nxcd = np.linspace(np.amin(nxgd),np.amax(nxgd),cfg.pms['shape'][0])
    phcd = np.linspace(np.amin(nygd),np.amax(nygd),cfg.pms['shape'][1])
    rcd = np.linspace(np.amin(rgd),np.amax(rgd),cfg.pms['shape'][2])
    ndotq_s = griddata((nxgd, nygd, rgd), ndotq, (nxgd, nygd, rcd), method='linear')
    ndotq=None;
    # rgd, nxgd, nygd, ndotq_s, tau
 # Compute the kSZ signal
    fn = np.exp(-tau)*ndotq_s/cfg.pms['c']*mToMpc*cfg.pms['Tcmb']*(1+z) # (Mpc/s)(m/s)^-1(m/Mpc)(K) = K
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
    dTkSZ = compute_kSZ2()
    plt.imshow(dTkSZ,origin='lower')
    plt.xlabel(r"$x\ (\mathrm{Mpc})$",fontsize=18)
    plt.ylabel(r"$y\ (\mathrm{Mpc})$",fontsize=18)
    cb=plt.colorbar()
    cb.set_label(r"$\Delta T_{kSZ}$",fontsize=18)
    plt.show()
    np.save('kSZ_array2',dTkSZ)

    # dTkSZ = np.load('kSZ_array.npy')
    # compute_kSZ_pspec(dTkSZ)

