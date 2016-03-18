import numpy as np
import os 
import fnmatch as fnm
import re
import pylab as plt
import matplotlib as m 
m = reload(m) # I was getting white lines in the imshow plots. this fixed it. idkw.
import seaborn as sns
import scipy as sp

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
    pms['H0']=pms['h']*3.241e-18 # s^-1
    pms['mp']=1.672622e-27 # kg  
    pms['G']=6.67384e-11 # N m^2 kg^-2
    pms['c']=2.99792458e8 # m/s
    pms['YBBNp']=0.25
    pms['mHe']=6.6464764e-27 # kg
    pms['mH']=1.6737236e-27 # kg
    pms['mp']=1.6726219e-27 # kg
    pms['mu']=1.0 + pms['YBBNp']*0.25*(pms['mHe']/pms['mH']-1.0)
    pms['nb0']=(3.*pms['H0']**2*pms['Omb'])/(8.*np.pi*pms['G']*pms['mu']*pms['mp']) # at z=0. scales with (1+z)^3
    pms['sigT']=6.65246e-29 # m^2
    return pms 

pms = get_params(run_dir)
pms = get_constants(pms)

# *** I think this is off by 2! ***
def redshift_to_space(zi=0, zf=20, num=1000):
    dz = (zf-zi)/num
    z = np.linspace(zi,zf,num=num)
    fn = (pms['c'] / pms['H0']) / np.sqrt(pms['Omm']*(1+z)**3+(1.-pms['Omm']))
    d = sp.integrate.cumtrapz(fn,z,dx=dz,initial=0) / mToMpc
    return z,d

def space_to_redshift(d,zi=5,zf=8):
    z0,d0 = redshift_to_space(zi=zi,zf=zf,num=10000)
    z = np.interp(d, d0, z0)
    return z

def get_box_data(field,redshift):
    data_file = data_dir+match_file(data_dir,'{0}z{1:06.2f}*Mpc'.format(box_headers[field],redshift))
    return np.fromfile(data_file,dtype=np.float32)

def plot_box_data(field,redshift):
    nf1 = get_box_data(field,redshift)
    chosenIndex=200
    nf1 = nf1.reshape((400,400,400))[chosenIndex,:,:]
    plt.imshow(nf1)
    plt.show()

def get_int_box_data(field):
    # check if concatenated box exists
    f = match_file(data_dir,'{0}*lighttravel_cat.npy'.format(box_headers[field]))
    if f != None:  
        zi,zf,_,box_sizez,box_sizexy = [float(i) for i in re.findall("[-+]?\d+[\.]?\d*",f)] # regex magic that extracts numbers from a string
        print zi, zf, box_sizez
        pms['zMpc'] = int(box_sizez); pms['xyMpc'] = int(box_sizexy)
        pms['zi'] = zi;         pms['zf'] = zf
        box = np.load(data_dir+f)
        return box
    # create and save concatenated box
    flist = match_files(data_dir,'{0}*lighttravel'.format(box_headers[field]))
    box_list = []; zi_list = []; zf_list = []; box_zMpc = 0;
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
    print 'catbox',catbox.shape
    np.save('{5}{0}zstart{1}_zend{2}_FLIPBOXES1_{3}_{4}Mpc_lighttravel_cat'.format(box_headers[field],pms['zi'],pms['zf'],int(box_zMpc),int(box_size),data_dir),catbox)
    return catbox

def plot_int_box_data(field):
    nf1 = get_int_box_data(field)
    chosenIndex=200
    print "nf1.shape = ",nf1.shape
    nf1_slice = nf1[chosenIndex,:,:]
    
    #cmap = sns.cubehelix_palette(8, start=.5, rot=-.75, as_cmap=True)
    #sns.color_palette("BrBG", 100,as_cmap=True)
    plt.imshow(nf1_slice,cmap='jet')
    #cb = plt.colorbar()
    #cb.set_label('neutral fraction')
    #plt.ylabel('y (Mpc)')
    #plt.xlabel('z (Mpc)')
    plt.show()

def get_density_weights():
    density = get_int_box_data('density')
    denWeights = np.zeros_like(density)
    for ii in range(density.shape[-1]):
        denWeights[:,:,ii] = density[:,:,ii] / np.mean(density[:,:,ii])
    return denWeights

def compute_tau(density,nf):
    z0,d0 = redshift_to_space(zi=pms['zi'],zf=pms['zf'],num=1000)
    d = np.linspace(d0[0],d0[-1],nf.shape[2])
    z = space_to_redshift(d,zi=pms['zi'],zf=pms['zf'])
    nb = pms['nb0']*(1+z)**3
    chiHII = np.average(nf*(1+density),axis=(0,1))
    plt.plot(z,chiHII)
    chiHeIII = np.zeros_like(chiHII)
    chiHeIII[np.where(z<3)] = 1
    fn = pms['sigT']*nb*(chiHII+0.25*chiHeIII*pms['YBBNp'])*mToMpc
    tau = sp.integrate.cumtrapz(fn,dx=(d[-1]-d[0])/len(d),initial=0)
    #dldz = (pms['c']*pms['H0'])/((1.+z)*np.sqrt(pms['Omm']*(1.+z)**3+(1.-pms['Omm'])))
    return z,d,tau


if __name__=='__main__':
    density = get_int_box_data('density')
    nf = get_int_box_data('nf')
    z,d,tau = compute_tau(density,nf)
    plt.plot(z,tau)
    plt.show()
    
 



