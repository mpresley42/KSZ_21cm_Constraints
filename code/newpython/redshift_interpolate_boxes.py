import sys
import numpy as np
import scipy.integrate
from magic import * 
#import pylab as plt

########################################
# Get Command-line Arguments
########################################

# parameter file (varyTvir.dat)
param_file = sys.argv[1] 

# get list of boxes from file
boxlist_fname = sys.argv[2]
flist = open(boxlist_fname,'r').read().splitlines()


########################################
# Get Cosmo Parameters
########################################

# open the parameter file
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
drcdz = lambda z: (pms['c']/pms['H0'])/np.sqrt(pms['Omm']*(1+z)**3+(1-pms['Omm']))
# dr/dz in proper coordinates
drpdz = lambda z: (pms['c']/pms['H0'])/((1+z)*np.sqrt(pms['Omm']*(1+z)**3+(1-pms['Omm'])))

def z_to_d(z, num=10000, proper=False, array=False):
    """Takes in a redshift and returns its corresponding 
    comoving (or proper) distances in Mpc."""
    dz = float(z)/num
    z_arr = np.linspace(0,z,num=num)
    if proper: fn = drpdz(z_arr)
    else:      fn = drcdz(z_arr)
    if array:
        d_arr = scipy.integrate.cumtrapz(fn,z_arr,dx=dz)/pms['mToMpc']
        return z_arr[:-1], d_arr
    else:
        d = np.trapz(fn,z_arr,dx=dz)/pms['mToMpc']
        return d

def d_to_z(d, zmax=100, num=10000, proper=False):
    """Takes in a comoving (or proper) distance in Mpc which should be before
    zmax and returns the corresponding redshift."""
    z0,d0 = z_to_d(zmax, num=num, proper=proper, array=True)
    z = np.interp(d, d0, z0)
    return z

########################################
# Initialize Box Distance/Redshift Data
########################################

# get a list of the redshifts of each box
flen = len(flist)
zlist = np.zeros((flen,))
for ii,boxfname in enumerate(flist):
    # get the redshift, box size, and box size in Mpc from the filename
    fname_params = nums_from_string(boxfname)
    zlist[ii] = float(fname_params[0])
    box_size = int(fname_params[-2])
    box_size_Mpc = int(fname_params[-1])
    box_shape = (box_size,box_size,box_size)
    box_shape_Mpc = (box_size_Mpc,box_size_Mpc,box_size_Mpc)
    prefix = boxfname.split('_z')[0]
print 'zlist = ',zlist

# get a list of the co-moving distances of each box
dclist = np.zeros((flen,box_size))
for ii,boxz in enumerate(zlist):
    dc0 = z_to_d(boxz)
    dcbox = dc0 + np.arange(-box_size_Mpc/2.,box_size_Mpc/2.,
        float(box_size_Mpc)/box_size)
    dclist[ii,:] = dcbox
print 'dclist[:,0] = ',dclist[:,0]

########################################
# Create the Interpolated Box
########################################

interp_box = np.zeros((box_size,box_size,flen*box_size))
interp_dc = np.linspace(dclist[0,0],dclist[-1,-1],flen*box_size)
print interp_dc[0],interp_dc[-1]
print (interp_dc[-1] - interp_dc[0])/(flen*box_size)

# plt.plot(interp_dc,np.ones_like(interp_dc),label='interp_dc')
# for ii in xrange(dclist.shape[0]):
#     plt.plot(dclist[ii,:],(ii+2)*np.ones_like(dclist[ii,:]),label='box {0}'.format(ii))
# plt.legend()
# plt.show()

# quit()

########################################
# Main Program Loop 
########################################

weight_sum = np.zeros_like(interp_dc)
interp_z = np.zeros_like(interp_dc)
# loop over all of the original files
for jj in xrange(flen):
    print 'Now on box ',jj
    # load the box
    boxj = np.fromfile(flist[jj],dtype=np.float32).reshape(box_shape)
    
    # loop over the interpolated box
    for kk in xrange(len(interp_dc)):
        
        # find the redshift at this slice of the interpolated box
        zk = d_to_z(interp_dc[kk])
        interp_z[kk] = zk

        # find index of closest slice, assuming the boxes are periodic
        close_id = kk%box_size #np.argmin(np.abs(dclist[jj,:]-interp_dc[kk]))

        # assign the weight
        weight = 1./np.abs(zk-zlist[jj]) #1./np.abs(interp_dc[kk] - dclist[jj,close_id])
        if np.isinf(weight): weight = 10.
        #print zk, weight

        # add slice to weighted average
        interp_box[:,:,kk] += weight*boxj[:,:,close_id]

        # add weight to the weight sum 
        weight_sum[kk] += weight

# divide by the sum of weights to get the weighted average
interp_box = interp_box / weight_sum

########################################
# Save Data 
########################################

# split interpolated data into chunks
interp_chunks = np.split(interp_box,flen,axis=2)
interp_box=None

# save each chunk separately
for ii,chunk in enumerate(interp_chunks):
    # create the filename for interpolated data
    # interp_fname = '{0}_zstart{1:.3f}_zend{2:.3f}_FLIPBOXES42_{3}_{4}_{5}_{6:.2f}_lighttravel'.format(
    #                 prefix,zlist[0],zlist[-1],interp_box.shape[0],interp_box.shape[2],
    #                 box_size_Mpc,interp_dc[-1]-interp_dc[0])
    interp_fname = '{0}_zstart{1:06.2f}_zend{2:06.2f}_FLIPBOXES42_{3}_{4}_{5:.2f}_lighttravel'.format(
                    prefix,interp_z[ii*box_size],interp_z[(ii+1)*box_size-1],box_size,
                    box_size_Mpc,interp_dc[(ii+1)*box_size-1]-interp_dc[ii*box_size])
    print 'Saved to file: ',interp_fname

    # save the interpolated box
    chunk.astype(np.float32).tofile(interp_fname)



