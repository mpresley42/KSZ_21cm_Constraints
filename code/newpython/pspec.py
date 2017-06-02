import numpy as np

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
    kmin = np.min(kg)#;kmax = np.max(kg)
    kmax = min((np.max(kx),np.max(ky)))
    #kmin = np.min(np.abs(kx));kmax = np.max(kx)
    dk = (kmax-kmin)/n
    kbins = np.arange(kmin,kmax+dk,dk)
    #kbins = np.linspace(0.0,10.0)
    pspec = np.zeros_like(kbins)
    area = np.zeros_like(kbins)
    for ii in range(len(kx)):
        for jj in range(len(ky)):
            kbini = int(np.floor(kg[ii,jj]/dk))
            if kbini < len(pspec):
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
    """Pad 2d array with zeros so that it extends to a certain size.
    Also extends coordinate arrays in a linear fashion to extend to 
    that size."""
    
    # if the pad width is less than the size of the array, then return
    # the original array
    if pdw <= len(nx): return nx,ny,nMap
    if pdw <= len(ny): return nx,ny,nMap
    
    # determine the amount to pad to the left and right
    # for the x direction
    pdx = (pdw - nMap.shape[0])/2.
    pdxl=int(np.ceil(pdx))
    pdxr=int(np.floor(pdx))

    # determine the amount to pad to the left and right 
    # for the y direction
    pdy = (pdw - nMap.shape[1])/2.
    pdyl=int(np.ceil(pdy))
    pdyr=int(np.floor(pdy))
    
    # pad the 2d array with zeros
    nMap_pad = np.pad(nMap, ((pdxl,pdxr),(pdxl,pdxr)), 'constant', constant_values=0)
    
    # pad the x coordinates with a linear extension
    dx=(nx[-1]-nx[0])/(len(nx)-1.); #print dx
    nx_pad = np.pad(nx,(pdxl,pdxr),'linear_ramp',end_values=(nx[0]-pdxl*dx,nx[-1]+(pdxr+1)*dx))
    
    # pad the y coordinates with a linear extension
    dy=(ny[-1]-ny[0])/(len(ny)-1.)
    ny_pad = np.pad(ny,(pdyl,pdyr),'linear_ramp',end_values=(ny[0]-pdyl*dy,ny[-1]+(pdyr+1)*dy))
    
    return nx_pad,ny_pad,nMap_pad
