import numpy as np
import scipy as sp
import scipy.integrate
import fnmatch as fnm
import re
import os 

box_headers = {'density':"updated_smoothed_deltax_",'vx':"updated_vx_",'vy':"updated_vy_",'vz':"updated_vz_",'nf':"xH_nohalos_"} 

def match_file(location,pattern):
    for file in os.listdir(location):
        if fnm.fnmatch(file, pattern):
            return file
    return None

def match_files(location,pattern):
    files = []
    for f in os.listdir(location):
        if fnm.fnmatch(f, pattern):
            files.append(f)
    return files

def nums_from_string(f):
    # regex magic that extracts numbers from a string
    return [float(i) for i in re.findall("[-+]?\d+[\.]?\d*",f)] 

#############################################
# One run of a 21cmFAST simulation
#############################################
class Sim21cm:
    def __init__(self,data_dir,run_dir,itot=8):
        self.data_dir = data_dir 
        self.run_dir = run_dir
        self.itot = itot
        self.pms = {}
        self.set_constants()
        self.get_params(run_dir)
        self.setup_iboxes()
        self.set_ibox_shape()

    ##########################
    # set up pms
    ##########################
    def get_params(self,run_dir):
     # return a dictionary of parameters based off of the directory name
        params=run_dir.split('_')
        for ii in xrange(0,len(params),2):
            self.pms[params[ii]] = float(params[ii+1])

    def set_constants(self):
        self.pms['Sigma8']=0.83
        self.pms['h']=0.67
        self.pms['Omm']=0.32
        self.pms['Omb']=0.022/(self.pms['h']**2)
        self.pms['YBBNp']=0.24
        # =============================
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

    def set_ibox_shape(self):
        f = match_file(self.data_dir,'{0}*lighttravel_cat_0_{1}.npy'.format(box_headers['nf'],self.itot))
        if f==None:
            raise RuntimeError('No file found for field nf, ibox 0 when getting ibox shape')
        else:  
            self.check_assign_box_pms(f)
            box = np.load(self.data_dir+f)
            self.pms['ishape'] = box.shape

    ##########################
    # set up box data
    ##########################
    def setup_iboxes(self):
        for field in box_headers.keys():
            # check if iboxes exist
            flist = match_files(self.data_dir,'{0}*lighttravel_cat_*_*.npy'.format(box_headers[field]))
            if flist == []:  
                self.create_iboxes(field)
            f = match_file(self.data_dir,'{0}*lighttravel_cat_0_*.npy'.format(box_headers[field]))            
            self.check_assign_box_pms(f)
            print "set up iboxes for field %s" % field
            
    def create_iboxes(self,field):
        # load all of the integrated lighttravel boxes
        flist = match_files(self.data_dir,'{0}*lighttravel'.format(box_headers[field]))
        box_list = []; zi_list = []; zf_list = []; box_zMpc = 0;
        if len(flist)==0: raise RuntimeError('No file found for field {0}'.format(field))
        for f in flist:
            zi,zf,_,box_size,box_Mpc = nums_from_string(f)
            data_file = self.data_dir+f
            box_data = np.fromfile(data_file,dtype=np.float32)
            box_list.append(box_data.reshape((box_size,box_size,box_size)))
            zi_list.append(zi)
            zf_list.append(zf)
            box_zMpc += box_Mpc
        # define parameters for concatenated box
        boxpms = {}
        boxpms['zMpc'] = int(box_zMpc)
        boxpms['xyMpc'] = int(box_Mpc)
        boxpms['zi'] = min(zi_list)
        boxpms['zf'] = max(zf_list)
        # create catbox
        sorted_box_list = [box for (zi,box) in sorted(zip(zi_list,box_list))]
        catbox = np.concatenate(sorted_box_list,axis=2)
        boxpms['shape']=catbox.shape
        # create and save separated iboxes 
        ibox_list = np.array_split(catbox,self.itot,axis=2)
        catbox=None
        for ii,box in enumerate(ibox_list):
            np.save('{5}{0}zstart{1}_zend{2}_FLIPBOXES1_{3}_{4}Mpc_lighttravel_cat_{6}_{7}'.format(
                box_headers[field],boxpms['zi'],boxpms['zf'],boxpms['zMpc'],boxpms['xyMpc'],self.data_dir,ii,self.itot),box)

    def parse_ibox_filename(self,f):
        zi,zf,_,box_zMpc,box_xyMpc,_,itot = nums_from_string(f)
        itot = int(itot)
        zMpc = int(box_zMpc)/itot; xyMpc = int(box_xyMpc)
        return zi,zf,zMpc,xyMpc,itot

    def check_assign_box_pms(self,f):
        boxpm_labels = ('zi','zf','zMpc','xyMpc','itot')
        boxpms = self.parse_ibox_filename(f)
        for lab,pm in zip(boxpm_labels,boxpms):
            if lab in self.pms: 
                assert self.pms[lab] == pm, "New {0} not consistent with loaded data.".format(lab)
            else:  
                self.pms[lab] = pm

    ##########################
    # redshift
    ##########################
    f = lambda self,z: (self.pms['c'] / self.pms['H0']) / np.sqrt(self.pms['Omm']*(1+z)**3+(1.-self.pms['Omm'])) # drdz co-moving
    fp = lambda self,z: (self.pms['c'] / self.pms['H0']) / ((1+z)*np.sqrt(self.pms['Omm']*(1+z)**3+(1.-self.pms['Omm']))) # drdz proper
    def _redshift_to_space(self,zi=0, zf=20, num=1000, proper=False):
        """Takes in a starting and ending redshift and 
           returns an array of num redshifts and an array 
           of their corresponding comoving distances in Mpc."""
        dz = (zf-zi)/num
        z0 = np.linspace(0,zi,num=num)
        if proper: fn0 = self.fp(z0)
        else:      fn0 = self.f(z0)    
        di = sp.integrate.trapz(fn0,z0,dx=dz)
        z = np.linspace(zi,zf,num=num)
        if proper: fn = self.fp(z)
        else:      fn = self.f(z)
        d = (di+sp.integrate.cumtrapz(fn,z,dx=dz,initial=0)) / self.pms['mToMpc']
        return z,d

    def _space_to_redshift(self,d,zi=5,zf=8,proper=False):
        """Takes in an array of comoving distances in Mpc and returns 
           the corresponding array of redshifts."""
        z0,d0 = self._redshift_to_space(zi=zi,zf=zf,num=10000,proper=proper)
        z = np.interp(d, d0, z0)
        return z

    def get_z_d(self,zi,zf,dlen=None,proper=False):
        """Takes in a starting and ending redshift and 
           returns arrays of redshifts and comoving distances  
           in Mpc for the coordinates of the boxes."""
        z0,d0 = self._redshift_to_space(zi=zi,zf=zf,num=1000,proper=proper)
        if dlen==None: d = np.linspace(d0[0],d0[-1],self.pms['ishape'][2])
        else: d = np.linspace(d0[0],d0[-1],dlen)
        z = self._space_to_redshift(d,zi=zi,zf=zf,proper=proper)
        return z,d

    def get_tot_z_d(self,proper=False):
        z,d = self.get_z_d(self.pms['zi'],self.pms['zf'],
                        dlen=self.pms['ishape'][2]*self.pms['itot'],proper=proper)
        return z,d    


#############################################
# Data boxes from a run of 21cmFAST 
#############################################
class Box21cm:
    def __init__(self,sim,ibox):
        self.sim = sim
        self.ibox = ibox
        self.itot = self.sim.itot
        self.ishape = self.sim.pms['ishape']

    ##########################
    # box loading
    ##########################
    def get_data(self,field):
        # check if ibox exists
        f = match_file(self.sim.data_dir,'{0}*lighttravel_cat_{1}_*.npy'.format(box_headers[field],self.ibox))
        if f==None:
            raise RuntimeError('No file found for field {0}, ibox {1}'.format(field,self.ibox))
        else:  
            self.sim.check_assign_box_pms(f)
            box = np.load(self.sim.data_dir+f)
            assert self.ishape == box.shape, "New ishape not consistent with loaded data."
            return box

    ##########################
    # box redshift
    ##########################
    def get_z_d(self,proper=False):
        z,d = self.sim.get_tot_z_d(proper=proper)
        zlist = np.array_split(z,self.itot)
        dlist = np.array_split(d,self.itot)
        return zlist[self.ibox],dlist[self.ibox]

    ##########################
    # box coordinates
    ##########################
    # Create a grid of x,y,z coordinates for the standard box
    def get_xyz_cd(self):
        box = self.sim.pms['ishape']
        boxMpc = np.array([self.sim.pms['xyMpc'],self.sim.pms['xyMpc'],self.sim.pms['zMpc']])
        xcd = np.arange(-boxMpc[0]/2,boxMpc[0]/2,boxMpc[0]/box[0])
        ycd = np.linspace(-boxMpc[1]/2,boxMpc[1]/2,boxMpc[1]/box[1])
        z,d = self.get_z_d(self.ibox)
        zcd = d[0] + np.arange(0,boxMpc[2],boxMpc[2]/box[2])
        return xcd,ycd,zcd
    def get_xyz_gd(self):
        xcd,ycd,zcd = self.get_xyz_cd()
        xyzgd = np.meshgrid(xcd,ycd,zcd)
        return xyzgd
    # Create 1D arrays of nx,ny,r coords for the standard box
    def get_nxnyr_cd(self):
        box = self.ishape # shape of ibox
        boxMpc = np.array([self.sim.pms['xyMpc'],self.sim.pms['xyMpc'],self.sim.pms['zMpc']]) # size of ibox in Mpc
        lx = boxMpc[0]/2.; ly = boxMpc[1]/2.; lz = boxMpc[2]*self.itot # size of total box in Mpc
        z,d = self.sim.get_tot_z_d() # z,d for total box
        zi,di = self.get_i_z_d() # z,d for ibox

        # back of box -- throws away half the box but whatever
        df = d[0]+lz # back of total box
        dif = di[0]+boxMpc[2] # back of ibox
        nx_max = lx / np.sqrt(lx*lx+df*df) # nx_min = - nx_max
        ny_max = ly / np.sqrt(ly*ly+df*df) # ny_min = - ny_max
        r_max = dif/np.sqrt(1. - nx_max*nx_max - ny_max*ny_max)
        r_min = di[0]/np.sqrt(1. - nx_max*nx_max - ny_max*ny_max) #np.sqrt(di[0]*di[0]+lx*lx+ly*ly) # r_min = d[0]

        nxcd = np.linspace(-nx_max,nx_max,box[0])
        nycd = np.linspace(-ny_max,ny_max,box[1])
        rcd = np.linspace(r_min,r_max,box[2])

        return nxcd,nycd,rcd




