import numpy as np
import h5py
import scipy as sp
import scipy.integrate
import pylab as plt
from magic import *
import timeit

HEADERS = {'density':"updated_smoothed_deltax_",'vx':"updated_vx_",'vy':"updated_vy_",'vz':"updated_vz_",'nf':"xH_nohalos_"} 

#HEADERS = {'density':"updated_smoothed_deltax_",}

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
        for ii in xrange(0,len(params),2):
            self.pms[params[ii]] = float(params[ii+1])

    def set_data_boxes(self):
        """Creates a dictionary of box objects that each wrap 
        the data in memmaps.
        https://docs.scipy.org/doc/numpy/reference/generated/numpy.memmap.html"""
        
        
        # Gets the parameters from the filenames of one of the fields
        # Note: all fields will have the same parameters
        zi_list = []; zf_list = []; box_zMpc = 0; num_boxes = 0
        flist = match_files(self.data_dir,'{0}*lighttravel'.format(HEADERS.items()[0][1]))
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

        # Create a base file name for the box object hdf5 files
        #updated_v_zstart005.00000_zend009.56801_FLIPBOXES0_1024_1600Mpc_lighttravel
        fbase = 'zstart{0}_zend{1}_FLIPBOXES0_{2}_{3}_{4}Mpc_{5}Mpc_lighttravel'.format(
            self.pms['zi'],self.pms['zf'],self.pms['shape'][0],
            self.pms['shape'][2],self.pms['xyMpc'],self.pms['zMpc'])
        
        # Create a dictionary containing the box objects for 
        # for all of the fields
        for field in HEADERS.keys():
            #Use magic.py to get the proper files for this box
            flist = match_files(self.data_dir,'{0}*lighttravel'.format(HEADERS[field]))
            fname = '{0}{1}'.format(HEADERS[field],fbase)
            self.box[field] = Box(self.data_dir,flist,fname,self.pms['shape'])


    ##########################
    # redshift
    ##########################
    
    # dr/dz in co-moving coordinates
    fc = lambda self,z: (self.pms['c'] / self.pms['H0']) / np.sqrt(self.pms['Omm']*(1+z)**3+(1.-self.pms['Omm']))
    # dr/dz in proper coordinates
    fp = lambda self,z: (self.pms['c'] / self.pms['H0']) / ((1+z)*np.sqrt(self.pms['Omm']*(1+z)**3+(1.-self.pms['Omm'])))
    
    def redshift_to_space(self,zi=0, zf=20, num=10000, proper=False):
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
    
    def space_to_redshift(self,d,zi=0,zf=20,proper=False):
        """Takes in an array of comoving distances in Mpc and returns 
           the corresponding array of redshifts."""
        z0,d0 = self.redshift_to_space(zi=zi,zf=zf,proper=proper)
        z = np.interp(d, d0, z0)
        return z

    def get_z_d(self,proper=False):
        """Gets arrays of redshift and comoving distances that 
           correspond to the coordinates of the box."""
        z0,d0 = self.redshift_to_space(zi=self.pms['zi'],
            zf=self.pms['zf'],proper=proper)
        d = np.linspace(d0[0],d0[-1],self.pms['shape'][2])
        z = self.space_to_redshift(d,zi=self.pms['zi'],
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

    # Create 1D arrays of nx,ny,r coords for the standard box
    def get_nxnyr_cd(self):
        """Get the angular nx, ny, r coordinates that correspond
        to the data box axes"""
        bsh = self.pms['shape'] 
        bMpc = np.array([self.pms['xyMpc'],self.pms['xyMpc'],
            self.pms['zMpc']])
        lx = bMpc[0]/2.; ly = bMpc[1]/2.; lz = bMpc[2]
        z,d = self.get_z_d() # z,d for total box

        # back of box -- throws away half the box but whatever
        df = d[0]+lz 
        nx_max = lx / np.sqrt(lx*lx+df*df) # nx_min = - nx_max
        ny_max = ly / np.sqrt(ly*ly+df*df) # ny_min = - ny_max
        r_max = df/np.sqrt(1. - nx_max*nx_max - ny_max*ny_max)
        r_min = d[0]/np.sqrt(1. - nx_max*nx_max - ny_max*ny_max) 

        nxcd = np.linspace(-nx_max,nx_max,bsh[0])
        nycd = np.linspace(-ny_max,ny_max,bsh[1])
        rcd = np.linspace(r_min,r_max,bsh[2])

        return nxcd,nycd,rcd


#############################################
# One Data Box (may span multiple files)
# This version uses h5py and is even slower than memmap
#############################################
class Box1:
    """This class hold h5py datasets containing the data box 
    that is stored in multiple files"""
    
    def __init__(self,data_dir,flist,fname,shape):
        """Creates the h5py dataset for the files"""
        self.data_dir = data_dir
        self.fname = fname

        h5file = h5py.File('{0}{1}'.format(data_dir,fname),"w")
        dset = h5file.create_dataset("data",shape,dtype='float32')
        arlist = []
        for ii,f in enumerate(flist):
            f0 = np.fromfile('{0}{1}'.format(data_dir,f), 
                dtype=np.float32).reshape(
                shape/np.array((1,1,len(flist))))
            starti = ii*shape[2]/len(flist)
            endi   = (ii+1)*shape[2]/len(flist)
            dset[:,:,starti:endi] = f0
            f0 = None
        h5file.close()

    def __getitem__(self,key):
        """Does the array thing, i.e. box[i,j,k]."""
        h5file = h5py.File('{0}{1}.hdf5'.format(self.data_dir,self.fname),"r")
        dset = h5file['data']
        return np.float32(dset[key])


    def slice(self,k):
        """Picks out one k-slice"""
        h5file = h5py.File('{0}{1}'.format(self.data_dir,self.fname),"r")
        dset = h5file['data']
        return dset[:,:,k]



#############################################
# One Data Box (may span multiple files)
# This version uses memmap and is very slow
#############################################
class Box2:
    """This class holds a memmap to a data box that is stored in
    multiple files."""

    def __init__(self,data_dir,flist,_,shape):
        """Creates memmaps to the files"""       
        # number of cubes inside box
        self.num_cubes = int(len(flist))
        # shape of box
        self.shape = shape
        self.cube_size = self.shape[2]/self.num_cubes
        # shape of cube
        self.cube_shape = (self.shape[0],self.shape[1],self.cube_size)

        self.cubes = []
        for f in flist:
            self.cubes.append(np.memmap('{0}{1}'.format(data_dir,f), 
                dtype=np.float32, mode='r', 
                shape=self.cube_shape))

    def __getitem__(self,key):
        """Does the array thing, i.e. box[i,j,k]. This will be slow,
        so it is recommended to use slice if possible."""
        icube = key[2]/self.cube_size
        ijkcube = (key[0],key[1],key[2]%self.cube_size)
        return self.cubes[icube][ijkcube]

    def slice(self,k):
        """Picks out one k-slice"""
        icube = int(k)/self.cube_size
        kcube = int(k)%self.cube_size
        return np.array(self.cubes[icube][:,:,kcube])

    def slicex(self,i):
        """Picks out one i-slice"""
        array_list = []
        for icube in xrange(self.num_cubes):
            array_list.append(np.array(self.cubes[icube][i,:,:]))
        return np.concatenate(array_list)

    def slicey(self,j):
        """Picks out one j-slice"""
        array_list = []
        for icube in xrange(self.num_cubes):
            array_list.append(np.array(self.cubes[icube][:,j,:]))
        return np.concatenate(array_list)


#############################################
# One Data Box (may span multiple files)
# This version chunks the files
#############################################
class Box:
    """This class splits the data box into chunks and only holds one
    chunk in memory at a time. If you ask for a part of the box not in 
    the current memory, it dumps the current chunk and replaces it with
    the new one."""

    def __init__(self,data_dir,flist,fname,shape,num_chunks=8):
        """Chunks the data and creates a new file
        for each chunk. Note: num_chunks must be divisible by
        the number of files """       
        # directory containing chunked files
        self.data_dir = data_dir
        # base file name for the chunk files
        self.fname = fname
        # list of file names for the original data boxes
        self.flist = flist
        # number of chunks inside box
        assert num_chunks%len(flist)==0
        self.num_chunks = num_chunks
        # shape of entire box
        self.shape = shape
        # shape of each chunk
        self.chunk_shape = (self.shape[0],self.shape[1],
            self.shape[2]/self.num_chunks)
        # chunk file name base
        self.fname = fname
        
        # create chunk files
        self._create_chunks()

        # place-holder for the currently loaded file
        self.current_chunk = None
        self.current_index = None

    def _create_chunks(self):
        """Create the chunked boxes and save the chunks in 
        their own files."""
        ii=0
        for f in self.flist:
            cube = np.fromfile('{0}{1}'.format(self.data_dir,f), 
                dtype=np.float32).reshape((self.shape[0],self.shape[1],
                    self.shape[2]/len(self.flist)))
            chunk_list = np.split(cube,self.num_chunks/len(self.flist),axis=2)
            for chunk in chunk_list:
                np.save('{0}{1}_{2}'.format(self.data_dir,self.fname,ii),chunk)
                ii+=1

    def __getitem__(self,key):
        """Does the array thing, i.e. box[i,j,k]. This will be slow,
        so it is recommended to use slice if possible."""
        ichunk = key[2]/self.chunk_shape[2]
        ijkchunk = (key[0],key[1],key[2]%self.chunk_shape[2])
        print 'index',ichunk, ijkchunk
        if ichunk!=self.current_index:
            self.current_index = ichunk
            self.current_chunk = np.load('{0}{1}_{2}.npy'.format(self.data_dir,self.fname,ichunk))
        return self.current_chunk[ijkchunk]


    def slice(self,k):
        ichunk = k/self.chunk_shape[2]
        kchunk = k%self.chunk_shape[2]
        if ichunk!=self.current_index:
            self.current_index = ichunk
            self.current_chunk = np.load('{0}{1}_{2}.npy'.format(self.data_dir,self.fname,ichunk))
        return self.current_chunk[:,:,kchunk]


if __name__=='__main__':
    #data_dir = '/Users/mpresley/Research/KSZ_21cm_Constraints/data/mesinger_1/original_1_mesinger_1/'
    data_dir = '/Users/mpresley/Research/KSZ_21cm_Constraints/data/small/RandSeed_111_Sigma8_0.81577_h_0.68123_Omm_0.30404_Omb_0.04805_ns_0.96670_Rmfp_35.00_Tvir_60000.0_Squiggly_40.00_lnAs_3.06400/'

    start = timeit.default_timer()
    sim = Sim(data_dir) 
    end = timeit.default_timer()
    print 'Simbox initialization: ', end-start

    sim.get_nxnyr_cd()

    # print sim.run_dir
    # print sim.pms['zi'] 
    # print sim.pms['zf'] 
    # print sim.pms['zMpc'] 
    # print sim.pms['xyMpc'] 
    # print sim.pms['shape']
    # rhobox = sim.box['density']
    # print rhobox[0,0,100]
    # print rhobox[0,0,200]
    # print rhobox.slice(200)[0,0]
    # print rhobox.slice(100).shape
    # z,d = sim.get_z_d()
    # plt.plot(d,z,c='orange')
    # plt.axhline(sim.pms['zi'],c='steelblue')
    # plt.axhline(sim.pms['zf'],c='forestgreen')
    # plt.show()
    # xyzcd = sim.get_xyz_coords()
    # print xyzcd[0].shape, xyzcd[1].shape, xyzcd[2].shape






