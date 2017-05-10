import numpy as np
import scipy as sp
import scipy.integrate
import pylab as plt
from magic import *

HEADERS = {'density':"updated_smoothed_deltax_",'vx':"updated_vx_",'vy':"updated_vy_",'vz':"updated_vz_",'nf':"xH_nohalos_"} 

#############################################
# One run of a 21cmFAST simulation
#############################################
class Sim:
    """This class holds all of the parameters and the data locations
       for a single run of 21cmFAST. It allows the user to smoothly """

    def __init__(self,data_dir):
        """Take in the location of the data and set all of the 
        simulation parameters and inialize the memmaps for the 
        data boxes."""
        self.data_dir = data_dir
        self.run_dir  = self.data_dir.split('/')[-2]
        self.pms = {}
        self.set_pms()
        self.box = {}
        self.set_data_boxes()

    def set_pms(self):
        """Sets cosmological constants and recovers simulation 
        parameters from name of the directory that 21cmFAST stores
        the data in. These are stored in the pms dictionary"""
        
        # Set the default parameters
        self.pms['Sigma8']=0.83
        self.pms['h']=0.67
        self.pms['Omm']=0.32
        self.pms['Omb']=0.022/(self.pms['h']**2)
        
        # Set the cosmological constants
        self.pms['YBBNp']=0.24
        self.pms['H0']=self.pms['h']*3.241e-18 # s^-1 # = h * 100 km/s/Mpc
        self.pms['mp']=1.672622e-27 # kg  
        self.pms['G']=6.67384e-11 # N m^2 kg^-2
        self.pms['c']=2.99792458e8 # m/s
        self.pms['mHe']=6.6464764e-27 # kg
        self.pms['mH']=1.6737236e-27 # kg
        self.pms['mp']=1.6726219e-27 # kg
        self.pms['mu']=1.0 + self.pms['YBBNp']*0.25*(
                             self.pms['mHe']/self.pms['mH']-1.0)
        self.pms['nb0']=(3.*self.pms['H0']**2*self.pms['Omb'])/(
                         8.*np.pi*self.pms['G']*self.pms['mu']*
                         self.pms['mp']) # at z=0. scales with (1+z)^3
        self.pms['sigT']=6.65246e-29 # m^2
        self.pms['Tcmb']=2.725e6 # microK
        self.pms['mToMpc'] = 3.086e22 # meters in a Mpc
        
        # Set the parameters from the directory name
        params=self.run_dir.split('_')
        print params
        for ii in xrange(0,len(params),2):
            self.pms[params[ii]] = float(params[ii+1])

    def set_data_boxes(self):
        """Creates a dictionary of box objects that each wrap 
        the data in memmaps.
        https://docs.scipy.org/doc/numpy/reference/generated/numpy.memmap.html"""
        
        for field in HEADERS.keys():
            #Use magic.py to get the proper files for this box
            flist = match_files(self.data_dir,'{0}*lighttravel'.format(HEADERS[field]))
            _,_,_,cube_size,_ = nums_from_string(flist[0])
            self.box[field] = Box(self.data_dir,flist,field,int(cube_size))
        
        # Gets the parameters from the filenames of the last field
        # Note: all fields will have the same parameters
        zi_list = []; zf_list = []; box_zMpc = 0; num_boxes = 0
        print flist
        for f in flist:
            zi,zf,_,box_size,box_Mpc = nums_from_string(f)
            zi_list.append(float(zi))
            zf_list.append(float(zf))
            box_zMpc += box_Mpc
        self.pms['zi'] = min(zi_list)
        self.pms['zf'] = max(zf_list)
        self.pms['zMpc'] = int(box_zMpc)
        self.pms['xyMpc'] = int(box_Mpc)
        self.pms['shape'] = (int(box_size),int(box_size),len(flist)*int(box_size))

    ##########################
    # redshift
    ##########################
    
    # dr/dz in co-moving coordinates
    fc = lambda self,z: (self.pms['c'] / self.pms['H0']) / np.sqrt(self.pms['Omm']*(1+z)**3+(1.-self.pms['Omm']))
    # dr/dz in proper coordinates
    fp = lambda self,z: (self.pms['c'] / self.pms['H0']) / ((1+z)*np.sqrt(self.pms['Omm']*(1+z)**3+(1.-self.pms['Omm'])))
    
    def _redshift_to_space(self,zi=0, zf=20, num=10000, proper=False):
        """Takes in a starting and ending redshift and 
           returns an array of num redshifts and an array 
           of their corresponding comoving distances in Mpc."""
        dz = (zf-zi)/num
        z0 = np.linspace(0,zi,num=num)
        if proper: fn0 = self.fp(z0)
        else:      fn0 = self.fc(z0)    
        di = sp.integrate.trapz(fn0,z0,dx=dz)
        z = np.linspace(zi,zf,num=num)
        if proper: fn = self.fp(z)
        else:      fn = self.fc(z)
        d = (di+sp.integrate.cumtrapz(fn,z,dx=dz,initial=0)) / self.pms['mToMpc']
        return z,d
    
    def _space_to_redshift(self,d,zi=5,zf=8,proper=False):
        """Takes in an array of comoving distances in Mpc and returns 
           the corresponding array of redshifts."""
        z0,d0 = self._redshift_to_space(zi=zi,zf=zf,proper=proper)
        z = np.interp(d, d0, z0)
        return z

    def get_z_d(self,proper=False):
        """Gets arrays of redshift and comoving distances that 
           correspond to the coordinates of the box."""
        z0,d0 = self._redshift_to_space(zi=self.pms['zi'],
            zf=self.pms['zf'],proper=proper)
        d = np.linspace(d0[0],d0[-1],self.pms['shape'][2])
        z = self._space_to_redshift(d,zi=self.pms['zi'],
            zf=self.pms['zf'],proper=proper)
        return z,d

    def get_xyz_coords(self):
        """Get the physical xyz distance coordinates that 
        correspond to the data box axes"""
        bsh = self.pms['shape']
        bMpc = np.array([self.pms['xyMpc'],self.pms['xyMpc'],
            self.pms['zMpc']])
        xcd = np.arange(-bMpc[0]/2.,bMpc[0]/2.,float(bMpc[0])/bsh[0])
        ycd = np.arange(-bMpc[1]/2.,bMpc[1]/2.,float(bMpc[1])/bsh[1])
        _,zcd = self.get_z_d()
        return xcd,ycd,zcd


#############################################
# One Data Box
#############################################
class Box:
    """This class holds a memmap to a data box that is stored in
    multiple files."""

    def __init__(self,data_dir,flist,field,cube_size):
        """Creates memmaps to the files"""       
        self.num_cubes = int(len(flist))
        self.cube_size = int(cube_size)
        self.cubes = []
        for f in flist:
            self.cubes.append(np.memmap('{0}{1}'.format(data_dir,f), 
                dtype=np.float32, mode='r', 
                shape=(self.cube_size,self.cube_size,self.cube_size)))

    def __getitem__(self,key):
        """Does the array thing, i.e. box[i,j,k]. This will be slow,
        so it is recommended to use slice if possible."""
        icube = key[2]/self.cube_size
        ijkcube = (key[0],key[1],key[2]%self.cube_size)
        print icube, ijkcube
        return self.cubes[icube][ijkcube]


    def slice(self,k):
        """Picks out one k-slice"""
        icube = int(k)/self.cube_size
        kcube = int(k)%self.cube_size
        return self.cubes[icube][:,:,kcube]

if __name__=='__main__':
    data_dir = '/Users/mpresley/Research/KSZ_21cm_Constraints/data/mesinger_1/original_1_mesinger_1/'
    sim = Sim(data_dir) 
    print sim.run_dir
    print sim.pms['zi'] 
    print sim.pms['zf'] 
    print sim.pms['zMpc'] 
    print sim.pms['xyMpc'] 
    print sim.pms['shape']
    rhobox = sim.box['density']
    print rhobox[0,0,1000]
    print rhobox[0,0,2000]
    print rhobox.slice(2000)[0,0]
    print rhobox.slice(1000).shape
    z,d = sim.get_z_d()
    # plt.plot(d,z,c='orange')
    # plt.axhline(sim.pms['zi'],c='steelblue')
    # plt.axhline(sim.pms['zf'],c='forestgreen')
    # plt.show()
    xyzcd = sim.get_xyz_coords()
    print xyzcd[0].shape, xyzcd[1].shape, xyzcd[2].shape






