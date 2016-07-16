import cfg
import numpy as np 
import scipy as sp
import scipy.integrate
from load_data import *
from redshift import *
from scipy.interpolate import griddata
from warnings import warn 

def compute_tau(density,nf,pretty=False):
    z0,d0 = get_z_d(0,cfg.pms['zi'],proper=True)
    z1,d1 = get_z_d(cfg.pms['zi'],cfg.pms['zf'],proper=True)
    z = np.concatenate((z0,z1))
    d = np.concatenate((d0,d1))

    nb = cfg.pms['nb0']*(1+z)**3
    chiHII = np.ones_like(z)
    chiHII[len(d0):] = np.average((1.-nf)*(1+density),axis=(0,1))
    chiHeIII = np.zeros_like(chiHII)
    chiHeIII[np.where(z<3)] = 1

    fn = cfg.pms['sigT']*nb*(chiHII+0.25*chiHeIII*cfg.pms['YBBNp'])*mToMpc
    tau = sp.integrate.cumtrapz(fn,x=d,initial=0)
    
    if pretty:
        plt.plot(z,chiHII,label='chiHII')
        plt.plot(z,chiHeIII,label='chiHeIII')
        plt.ylim([0,1.1])
        plt.xlabel("z"); plt.ylabel("chi")
        plt.legend(loc='lower left'); plt.show() #bbox_to_anchor=(1.3,0.6)
        
        plt.plot(z,tau)
        plt.xlabel('z'); plt.ylabel('tau'); plt.show()

    # NOTE: This returns the proper distance!
    return z1,d1,tau[len(d0):]

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

# Create a grid of x,y,z coordinates for the standard box
def get_xyz_gd():
    box = cfg.pms['shape']
    xcd = np.linspace(-box[0]/2,box[0]/2,num=box[0])
    ycd = np.linspace(-box[1]/2,box[1]/2,num=box[1])
    z,d = get_z_d(cfg.pms['zi'],cfg.pms['zf'])
    zcd = d[0] + np.linspace(-box[2]/2,box[2]/2,num=box[2])
    xyzgd = np.meshgrid(xcd,ycd,zcd)
    xcd,ycd,zcd = None,None,None
    return xyzgd

# Create 1D arrays of nx,ny,r coords for the standard box
def get_nxnyr_cd():
    box = cfg.pms['shape']
    lx = box[0]/2.; ly = box[1]/2.; lz = box[2]
    z,d = get_z_d(cfg.pms['zi'],cfg.pms['zf'])
    
    # front of box -- don't use bc grid will extend
    #                 outside the standard box
    # nx_max = lx / np.sqrt(lx*lx+d[0]*d[0]) # nx_min = - nx_max
    # ny_max = ly / np.sqrt(ly*ly+d[0]*d[0]) # ny_min = - ny_max
    # r_max = np.sqrt(lx*lx+ly*ly+(d[0]+lz)*(d[0]+lz)) # r_min = d[0]

    # back of box -- throws away half the box but whatever
    df = d[0]+lz
    nx_max = lx / np.sqrt(lx*lx+df*df) # nx_min = - nx_max
    ny_max = ly / np.sqrt(ly*ly+df*df) # ny_min = - ny_max
    r_max = np.sqrt(lx*lx+ly*ly+df*df) # r_min = d[0]

    print nx_max, ny_max

    nxcd = np.linspace(-nx_max,nx_max,box[0])
    nycd = np.linspace(-ny_max,ny_max,box[1])
    print 2*nx_max/box[0], 2*ny_max/box[1]
    rcd = np.linspace(d[0],r_max,box[2])
    return nxcd,nycd,rcd

# Create a regular grid of nx,ny,r coords for the standard box
def get_nxnyr_gd():
    nxcd,nycd,rcd = get_nxnyr_cd()
    ngd = np.meshgrid(nxcd,nycd,rcd)
    return ngd

# Find the xyz values that correspond to the nx,ny,r coordinate grid
def get_corresp_xyz_gd():
 # create new nx,ny,r coordinate grids
    ngd = get_nxnyr_gd()
 # find the xyz values that correspond to the nx,ny,r coordinate grid
    n_xgd = ngd[0]*ngd[2]
    n_ygd = ngd[1]*ngd[2]
    n_zgd = np.sqrt(ngd[2]*ngd[2] - n_xgd*n_xgd - n_ygd*n_ygd)
    ngd=None
    return n_xgd,n_ygd,n_zgd

# Use nearest-neighbor to create a box on an nx,ny,nz grid
def regrid(ndotq):
 # find the xyz values that correspond to the nx,ny,r coordinate grid
    n_xgd,n_ygd,n_zgd = get_corresp_xyz_gd()
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
    n_igd=None; n_jgd=None; n_kgd=None
    if True: np.save('ndotq_ncd',ndotq_ncd)
    return ndotq_ncd

# Use linear interpolation to create a box on an nx,ny,nz grid
def regrid_linear(ndotq):
 # find the xyz values that correspond to the nx,ny,r coordinate grid
    n_xgd,n_ygd,n_zgd = get_corresp_xyz_gd()
 # find the indices of the closest (floor) point on the regular xyz grid
    z,d = get_z_d(cfg.pms['zi'],cfg.pms['zf'])
    box = np.array(cfg.pms['shape'])
    xyz_0 = -box/2.+np.array([0,0,d[0]])
    xyz_n = box/2.+np.array([0,0,d[0]])
    xyz_delta = (xyz_n - xyz_0)/box
    n_igd = ((n_xgd - xyz_0[0])/xyz_delta[0]).floor().clip(0,box[0]-1).astype(int)
    n_jgd = ((n_ygd - xyz_0[1])/xyz_delta[1]).floor().clip(0,box[1]-1).astype(int)
    n_kgd = ((n_zgd - xyz_0[2])/xyz_delta[2]).floor().clip(0,box[2]-1).astype(int)
 # create the interpolated grid of ndotq
 # y_* = y_j + (y_j+1 - y_j)(x_* - x_j)/(x_j+1 - x_j)
    # ndotq_ncd = ndotq[n_igd,n_jgd,n_kgd] + 
    #             (ndotq[n_igd+1,n_jgd+1,n_kgd+1]-ndotq[n_igd,n_jgd,n_kgd])*
    #             (n_xgd)

def compute_kSZ_linear(ndotq=None):
    if ndotq==None: ndotq = compute_ndotq()
    print "Have ndotq!"
    # ndotq    
 # Get tau
    nf = get_int_box_data('nf')
    density = get_int_box_data('density')
    zred,dp,tau = compute_tau(density,nf)
    zred,d = get_z_d(cfg.pms['zi'],cfg.pms['zf'],proper=False)
    # ndotq, nf, density, tau
    nf=None; density=None
    # ndotq, tau
    # Get n,r coord
    nxcd,nycd,rcd = get_nxnyr_cd()
    # Loop over angles
    dTkSZ = np.zeros([len(nxcd),len(nycd)])
    for ii in range(len(nxcd)):
        print ii
        for jj in range(len(nycd)):
            # Loop over z direction
            for kk in range(len(rcd)-1):
                # Find corresponding xyz coord
                x = nxcd[ii]*rcd[kk]
                y = nycd[jj]*rcd[kk]
                z = np.sqrt(rcd[kk]*rcd[kk]-x*x-y*y)
                # Find closest xyz indices
                dx=1.;dy=1;dz=1
                xi = int(np.floor(x/dx))
                yj = int(np.floor(y/dy))
                zk = int(np.floor((z-d[0])/dz))
                try:
                    dTkSZ[ii,jj] += np.exp(-tau[zk])*ndotq[xi,yj,zk]*(1+zred[zk])*(tau[zk+1]-tau[zk])
                except IndexError:
                    # print "Index Error"
                    # print ii,jj,kk
                    # print x,y,z 
                    # print xi,yj,zk
                    pass
    dTkSZ = dTkSZ/cfg.pms['c']*mToMpc*cfg.pms['Tcmb']
    return dTkSZ 

def compute_kSZ_test(ndotq=None):
    if ndotq==None: ndotq = compute_ndotq()
    print "Have ndotq!"
    # ndotq    
 # Get tau
    nf = get_int_box_data('nf')
    density = get_int_box_data('density')
    zred,d,tau = compute_tau(density,nf)
    # ndotq, nf, density, tau
    nf=None; density=None
    # ndotq, tau
    # z coordinate axis
    box = cfg.pms['shape']
    zcd = d[0] + np.linspace(-box[2]/2,box[2]/2,num=box[2])
    # Loop over angles
    dTkSZ = 0
    ii=200; jj=200
    # Loop over z direction
    for kk in range(box[2]-1):
        try:
            dTkSZ += np.exp(-tau[kk])*ndotq[ii,jj,kk]*(1+zred[kk])*(tau[kk+1]-tau[kk])
        except IndexError:
            print "Index Error: ",kk
            pass
    dTkSZ = dTkSZ/cfg.pms['c']*mToMpc*cfg.pms['Tcmb']
    return dTkSZ

def compute_kSZ(ndotq=None):
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
def pspec_2d(kx,ky,ft,n=500):
    ft = np.abs(ft)**2
    kxg,kyg = np.meshgrid(kx,ky)
    kg = np.sqrt(kxg*kxg+kyg*kyg)
    kmin = np.min(kg);kmax = np.max(kg)
    #kmin = np.min(np.abs(kx));kmax = np.max(kx)
    dk = (kmax-kmin)/n
    kbins = np.arange(kmin,kmax+dk,dk)
    #kbins = np.linspace(0.0,10.0)
    pspec = np.zeros_like(kbins)
    area = np.zeros_like(kbins)
    for ii in range(len(kx)):
        for jj in range(len(ky)):
            kbini = int(np.floor(kg[ii,jj]/dk))
            pspec[kbini] += ft[ii,jj]
            area[kbini] += 1
    pspec = pspec / area
    kbins=kbins[:-1]; pspec=pspec[:-1]; area=area[:-1]
    return kbins,pspec,area

def pad_array(nx,ny,nMap,pdw):
    nMap_pad = np.pad(nMap, ((pdw,pdw),(pdw,pdw)), 'constant', constant_values=0)
    dx=(nx[-1]-nx[0])/(len(nx)-1)
    nx_pad = np.concatenate((np.arange(nx[0]-pdw*dx,nx[0],dx),nx,np.arange(nx[-1],nx[-1]+pdw*dx,dx)))
    dy=(ny[-1]-ny[0])/(len(ny)-1)
    ny_pad = np.concatenate((np.arange(ny[0]-pdw*dy,ny[0],dy),ny,np.arange(ny[-1],ny[-1]+pdw*dy,dy)))
    return nx_pad,ny_pad,nMap_pad

def pad2_array(nx,ny,nMap,pdw):
    if pdw <= len(nx): return nx,ny,nMap
    if pdw <= len(ny): return nx,ny,nMap
    pdx = (pdw - nMap.shape[0])/2.; pdxl=np.ceil(pdx); pdxr=np.floor(pdx)
    pdy = (pdw - nMap.shape[1])/2.; pdyl=np.ceil(pdy); pdyr=np.floor(pdy)
    nMap_pad = np.pad(nMap, ((pdxl,pdxr),(pdxl,pdxr)), 'constant', constant_values=0)
    dx=(nx[-1]-nx[0])/(len(nx)-1.); print dx
    #nx_pad = np.concatenate((np.arange(nx[0]-pdxl*dx,nx[0],dx),nx,np.arange(nx[-1]+dx,nx[-1]+(pdxr+1)*dx,dx)))
    nx_pad = np.pad(nx,(pdxl,pdxr),'linear_ramp',end_values=(nx[0]-pdxl*dx,nx[-1]+(pdxr+1)*dx))
    dy=(ny[-1]-ny[0])/(len(ny)-1.)
    #ny_pad = np.concatenate((np.arange(ny[0]-pdyl*dy,ny[0],dy),ny,np.arange(ny[-1]+dy,ny[-1]+(pdyr+1)*dy,dy)))
    ny_pad = np.pad(ny,(pdyl,pdyr),'linear_ramp',end_values=(ny[0]-pdyl*dy,ny[-1]+(pdyr+1)*dy))
    return nx_pad,ny_pad,nMap_pad

def compute_kSZ_pspec(dTkSZ,mask=None,pdw=512,n=500,pretty=True):
 # Get the nx, ny, r coords
    nxcd,nycd,rcd = get_nxnyr_cd()
 # Apply mask if necessary 
    if mask==None: mask = np.ones_like(dTkSZ)
    else: dTkSZ = dTkSZ*mask
 # Compute overall normalization
    #print mask.shape, nycd.shape,nxcd.shape, dTkSZ.shape
    norm = np.trapz(np.trapz(mask*mask,nycd,axis=1),nxcd,axis=0)
    #print "norm = ",norm
 # Apply padding if necessary
    #if pdw!=0: nxcd,nycd,dTkSZ = pad_array(nxcd,nycd,dTkSZ,pdw)
    nxcd,nycd,dTkSZ = pad2_array(nxcd,nycd,dTkSZ,pdw)
 # Plot the new dTkSZ 
    if pretty:
        plt.imshow(dTkSZ,origin='lower')
        cb=plt.colorbar();cb.set_label(r"$\Delta T_{kSZ}$",fontsize=18)
        plt.show()
 # Find the Fourier Transform
    lx,ly,dTkSZ_FT = fft_2d(nxcd,nycd,dTkSZ) # [K][Mpc]^2    
    if pretty:
        plt.imshow(np.real(dTkSZ_FT),origin='lower')
        cb=plt.colorbar();cb.set_label(r"$Real(\Delta \widetildeT_{kSZ})$",fontsize=18)
        plt.show()
 # Find the Power Spectrum
    lbins,dTkSZ_P,area = pspec_2d(lx,ly,dTkSZ_FT,n=n) # [K]^2[Mpc]^4
    dTkSZ_P = dTkSZ_P/norm
    if pretty:
        plt.semilogy(lbins,dTkSZ_P)
        plt.ylabel(r"$P(\ell)$",fontsize=18)
        plt.xlabel(r"$\ell$",fontsize=18)
        plt.show()
    if pretty:
        plt.semilogy(lbins,lbins*(lbins+1.)*dTkSZ_P/(2*np.pi)) # [Mpc]^-2[K]^2[Mpc]^4 = [K]^2[Mpc]^2
        plt.ylabel(r"$\Delta_{kSZ}^2=\ell(\ell+1)C_\ell/2\pi [\mu K^2]$",fontsize=18)
        plt.xlabel(r"$k$",fontsize=18)
        plt.show()
    return lbins,dTkSZ_P,area

if __name__=='__main__':
    # test for compute_tau()
    
    nf = get_int_box_data('nf')
    
    print cfg.pms['zi'], cfg.pms['zf']

    density = get_int_box_data('density')
    #z,d,tau = compute_tau(density,nf,pretty=False)
    #plt.plot(z,tau); plt.show()
    #for kk in range(nf.shape[2]): print nf[200,200,kk]
    #plt.plot(z,nf[200,200]); plt.show()

    #get_nxnyr_cd()

    # ndotq = compute_ndotq()
    # np.save('ndotq_propercat',ndotq)
    # print "Got ndotq!"

    #ndotq = np.load('ndotq.npy')
    #dTkSZ = compute_kSZ_linear(ndotq)
    # print "dTkSZ = ",dTkSZ
    #np.save('dTkSZ_propercat',dTkSZ)
    
    dTkSZ = np.load('dTkSZ_propercat.npy')
    print np.mean(np.abs(dTkSZ))
    
    #print dTkSZ.shape

    plt.imshow(dTkSZ,origin='lower')
    plt.xlabel(r"$x\ (\mathrm{Mpc})$",fontsize=18)
    plt.ylabel(r"$y\ (\mathrm{Mpc})$",fontsize=18)
    cb=plt.colorbar()
    cb.set_label(r"$\Delta T_{kSZ}$",fontsize=18)
    plt.show()
    plt.hist(dTkSZ)
    plt.show()

    #compute_kSZ_pspec(dTkSZ,pdw=100,pretty=True)

    # LINEAR INTERPOLATION

    # dTkSZ = np.load('kSZ_array.npy')
    # compute_kSZ_pspec(dTkSZ)

    # also look at cube rotating thing in 21cm fast
    # ask Nick to save a sample of the cubes he is running
    # compare to reionization history
    # look at Messenger paper and look at the parameters he used
