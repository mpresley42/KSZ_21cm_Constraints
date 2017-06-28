import sys
import numpy as np
import scipy as sp
from magic import * 

########################################
# Parameters
########################################

# open the parameter file (varyTvir.dat)
param_file = sys.argv[1]
params = open(param_file,'r').readline()
#RandSeed Sigma8 h Omm Omb ns Rmfp Tvir Squiggly lnAs
_,_,h,Omm,_,_,_,_,_,_ = nums_from_string(params)

# set cosmological parameters
pms = {}
pms['c']      = 2.99792458e8 # m/s
pms['h']      = h
pms['H0']     = pms['h']*3.241e-18 # s^-1 # = h * 100 km/s/Mpc
pms['Omm']    = Omm
pms['mToMpc'] = 3.086e22 # meters in a Mpc

########################################
# Redshift-Space Conversion
########################################

# dr/dz in co-moving coordinates
dRdz = lambda z: (pms['c']/pms['H0'])/np.sqrt(pms['Omm']*(1+z)**3+(1-pms['Omm']))

def z_to_d(z, num=10000):
    """Takes in a redshift and returns its corresponding 
    comoving distances in Mpc."""
    dz = float(zf)/num
    z_arr = np.linspace(0,z,num=num)
    fn = dRdz(z_arr)
    d = sp.integrate.trapz(fn,z_arr,dx=dz,initial=0)/pms['mToMpc']
    return d

def d_to_z(d, zmax=100, num=10000):
    """Takes in a comoving distance in Mpc which should be before
    zmax and returns the corresponding redshift."""
    z0,d0 = z_to_d(zmax, num=num)
    z = np.interp(d, d0, z0)
    return z

########################################
# Initialize First Box
########################################

# open the box list file
boxlist_fname = sys.argv[2]
flist = open(boxlist_fname,'r').read().splitlines()

# get the first box's name
boxfname = flist[0]

# get the prefix of the filename
prefix = boxfname.split('_z')[0]

# get the redshift, box size, and box size in Mpc from the filename
fname_params = nums_from_string(boxfname)
z1 = float(fname_params[0])
box_size = int(fname_params[-2])
box_size_Mpc = int(fname_params[-1])
box_shape = (box_size,box_size,box_size)
box_shape_Mpc = (box_size_Mpc,box_size_Mpc,box_size_Mpc)

# read in the first box
box1 = np.fromfile(boxfname,dtype=np.float32).reshape(box_shape)

# initialize the interpolated box
#interp_box_shape = (box_size,box_size,box_size*len(flist))
#interp_box_shape_Mpc = (box_size_Mpc,box_size_Mpc,box_size_Mpc*len(flist))
interp_box_shape = box_shape
interp_box_shape_Mpc = box_shape_Mpc
interp_box = np.zeros(interp_box_shape,dtype=np.float32)
print 'interp_box_shape = ',interp_box_shape
print 'interp_box_shape_Mpc = ',interp_box_shape_Mpc


# set the spatial increment of all the boxes
dR = float(box_size_Mpc)*pms['mToMpc']/box_size # meters
print 'dR = {0} Mpc'.format(float(box_size_Mpc)/box_size)

# set the initial z of the interpolated box
zstart = z1
z = zstart

########################################
# Main Program Loop 
########################################

kk = 0
# loop over all of the boxes
for boxfname in flist[1:]:
    # get the second redshift
    z2 = nums_from_string(boxfname)[0]
    print "z1 = {0}, z2 = {1}".format(z1,z2)
    
    # read in the second box
    box2 = np.fromfile(boxfname,dtype=np.float32).reshape(box_shape)

    # increment over z until next box pair
    while z<z2:
        # increment z
        z += dR/dRdz(z) 
        # print 'dRdz = ',dRdz(z)
        # print "z = ",z

        # find closest slice from box1
        s1 = box1[:,:,kk]

        # find closest slice from box2
        s2 = box2[:,:,kk]

        # add linear interpolation between slices to the interpolated box
        interp_box[:,:,kk] = (s2-s1)*((z-z1)/(z2-z1))+s1
        #print interp_box[100,100,kk]
        kk += 1
        # check if we've filled up the interpolated box
        if kk==interp_box_shape[-1]:
            # write out the box 
            zend = z
            interp_fname = '{0}_zstart{1:.3f}_zend{2:.3f}_FLIPBOXES42_{3}_{4}_lighttravel'.format(
                prefix,zstart,zend,box_size,box_size_Mpc)
            interp_box.tofile(interp_fname)
            print "Wrote out interpolated box with zstart = {0:.3f} and zend = {1:.3f}".format(zstart,zend)
            interp_box = np.zeros(interp_box_shape,dtype=np.float32)
            # reset counts
            zstart = zend 
            kk = 0
    # make the second box the first box
    z1 = z2
    box1 = box2



