import cfg
import numpy as np 
import scipy as sp
import scipy.integrate
from load_data import *

f = lambda z: (cfg.pms['c'] / cfg.pms['H0']) / np.sqrt(cfg.pms['Omm']*(1+z)**3+(1.-cfg.pms['Omm'])) # drdz co-moving
fp = lambda z: (cfg.pms['c'] / cfg.pms['H0']) / ((1+z)*np.sqrt(cfg.pms['Omm']*(1+z)**3+(1.-cfg.pms['Omm']))) # drdz proper
def redshift_to_space(zi=0, zf=20, num=1000, proper=False):
    """Takes in a starting and ending redshift and 
       returns an array of num redshifts and an array 
       of their corresponding comoving distances in Mpc."""
    dz = (zf-zi)/num
    z0 = np.linspace(0,zi,num=num)
    if proper: fn0 = fp(z0)
    else:      fn0 = f(z0)    
    di = sp.integrate.trapz(fn0,z0,dx=dz)
    z = np.linspace(zi,zf,num=num)
    if proper: fn = fp(z)
    else:      fn = f(z)
    d = (di+sp.integrate.cumtrapz(fn,z,dx=dz,initial=0)) / mToMpc    
    return z,d

def space_to_redshift(d,zi=5,zf=8,proper=False):
    """Takes in an array of comoving distances in Mpc and returns 
       the corresponding array of redshifts."""
    z0,d0 = redshift_to_space(zi=zi,zf=zf,num=10000,proper=proper)
    z = np.interp(d, d0, z0)
    return z

def get_z_d(zi,zf,dlen=None,proper=False):
    """Takes in a starting and ending redshift and 
       returns arrays of redshifts and comoving distances  
       in Mpc for the coordinates of the boxes."""
    z0,d0 = redshift_to_space(zi=zi,zf=zf,num=1000,proper=proper)
    if dlen==None: d = np.linspace(d0[0],d0[-1],cfg.pms['shape'][2])
    else: d = np.linspace(d0[0],d0[-1],dlen)
    z = space_to_redshift(d,zi=zi,zf=zf,proper=proper)
    return z,d

if __name__=='__main__':
  nf = get_int_box_data('nf')
  print cfg.pms['shape']
  z,d = get_z_d(5.5,30.0,dlen=6722)
  zp = np.loadtxt('/Users/mpresley/Research/KSZ_21cm_Constraints/code/21cmFast/redshift_interpolate_z_log')
  plt.plot(d,z-zp)
  plt.show()

  z,d = get_z_d(0,10); plt.plot(z,d); 
  plt.plot((5.5,5.5),(0,10000)); plt.plot((8.3,8.3),(0,10000)) 
  plt.grid(b=True, which='major', color='k')
  plt.show()
