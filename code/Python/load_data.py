import numpy as np
import os 
import fnmatch as fnm
import re
# import scipy as sp
# import scipy.integrate
# import scipy.interpolate
import pylab as plt

run_dir='RandSeed_111_Sigma8_0.81577_h_0.68123_Omm_0.30404_Omb_0.04805_ns_0.96670_Rmfp_35.00_Tvir_60000.0_Squiggly_40.00_lnAs_3.06400'
data_dir='/Users/mpresley/Research/KSZ_21cm_Constraints/data/'+run_dir+'/'

box_headers = {'density':"updated_smoothed_deltax_",'vx':"updated_vx_",'vy':"updated_vy_",'vz':"updated_vz_",'v':"updated_v_",'nf':"xH_nohalos_"} 

mToMpc = 3.086e22 # meters in a Mpc
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
# return a dictionary of parameters based off of the directory name
def get_params(run_dir):
    pms=run_dir.split('_')
    params={}; ii=0
    while ii < len(pms):
        params[pms[ii]] = float(pms[ii+1])
        ii+=2
    return params

def get_constants(params):
    pms=params.copy()
    pms['Sigma8']=0.83
    pms['h']=0.67
    pms['Omm']=0.32
    pms['Omb']=0.022/(pms['h']*pms['h'])
    pms['YBBNp']=0.24
    # =============================
    pms['H0']=pms['h']*3.241e-18 # s^-1
    pms['mp']=1.672622e-27 # kg  
    pms['G']=6.67384e-11 # N m^2 kg^-2
    pms['c']=2.99792458e8 # m/s
    #pms['YBBNp']=0.25
    pms['mHe']=6.6464764e-27 # kg
    pms['mH']=1.6737236e-27 # kg
    pms['mp']=1.6726219e-27 # kg
    pms['mu']=1.0 + pms['YBBNp']*0.25*(pms['mHe']/pms['mH']-1.0)
    pms['nb0']=(3.*pms['H0']**2*pms['Omb'])/(8.*np.pi*pms['G']*pms['mu']*pms['mp']) # at z=0. scales with (1+z)^3
    pms['sigT']=6.65246e-29 # m^2
    pms['Tcmb']=2.725 # K
    return pms 

def make_pms():
    pms = get_params(run_dir)   
    pms = get_constants(pms)
    return pms

def get_box_data(field,redshift):
    data_file = data_dir+match_file(data_dir,'{0}z{1:06.2f}*Mpc'.format(box_headers[field],redshift))
    pms = make_pms()
    return np.fromfile(data_file,dtype=np.float32),pms

def plot_box_data(field,redshift):
    nf1, pms1 = get_box_data(field,redshift)
    chosenIndex=200
    nf1 = nf1.reshape((400,400,400))[chosenIndex,:,:]
    plt.imshow(nf1)
    plt.show()

def get_int_box_data(field):
    pms = make_pms()
    # check if concatenated box exists
    f = match_file(data_dir,'{0}*lighttravel_cat.npy'.format(box_headers[field]))
    if f != None:  
        zi,zf,_,box_sizez,box_sizexy = [float(i) for i in re.findall("[-+]?\d+[\.]?\d*",f)] # regex magic that extracts numbers from a string
        print zi, zf, box_sizez
        pms['zMpc'] = int(box_sizez); pms['xyMpc'] = int(box_sizexy)
        pms['zi'] = zi; pms['zf'] = zf
        box = np.load(data_dir+f)
        pms['shape']=box.shape
        return box, pms
    # create and save concatenated box
    flist = match_files(data_dir,'{0}*lighttravel'.format(box_headers[field]))
    box_list = []; zi_list = []; zf_list = []; box_zMpc = 0;
    if len(flist)==0: raise RuntimeError('No file found for field {0}'.format(field))
    for f in flist:
        print f
        zi,zf,_,box_size,_ = [float(i) for i in re.findall("[-+]?\d+[\.]?\d*",f)] # regex magic that extracts numbers from a string
        print zi, zf, box_size
        data_file = data_dir+f
        box_data = np.fromfile(data_file,dtype=np.float32)
        print "box_data.shape = ",box_data.shape
        #box_size = int(f.split('_')[-3])
        box_list.append(box_data.reshape((box_size,box_size,box_size)))
        zi_list.append(zi)
        zf_list.append(zf)
        box_zMpc += box_size
    pms['zMpc'] = int(box_zMpc)
    pms['xyMpc'] = int(box_size)
    pms['zi'] = min(zi_list)
    pms['zf'] = max(zf_list)
    sorted_box_list = [box for (zi,box) in sorted(zip(zi_list,box_list))]
    catbox = np.concatenate(sorted_box_list,axis=2)
    pms['shape']=catbox.shape
    print 'catbox',catbox.shape
    np.save('{5}{0}zstart{1}_zend{2}_FLIPBOXES1_{3}_{4}Mpc_lighttravel_cat'.format(box_headers[field],pms['zi'],pms['zf'],int(box_zMpc),int(box_size),data_dir),catbox)
    return catbox, pms 

def plot_int_box_data(field):
    nf1, pms1 = get_int_box_data(field)
    chosenIndex=200
    print "nf1.shape = ",nf1.shape
    nf1_slice = nf1[chosenIndex,:,:]
    plt.imshow(nf1_slice,cmap='jet')
    plt.show()

if __name__=='__main__':
  # nf,pms = get_int_box_data('nf')
  plot_box_data('nf',6.0)
  #plot_int_box_data('nf')
