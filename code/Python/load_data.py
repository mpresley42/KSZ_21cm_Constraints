import cfg
import numpy as np
import os 
import fnmatch as fnm
import re
import pylab as plt

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

def get_box_data(field,redshift):
    data_file = cfg.data_dir+match_file(cfg.data_dir,'{0}z{1:06.2f}*Mpc'.format(cfg.box_headers[field],redshift))
    return np.fromfile(data_file,dtype=np.float32)

def plot_box_data(field,redshift):
    nf1 = get_box_data(field,redshift)
    chosenIndex=200
    nf1 = nf1.reshape((400,400,400))[chosenIndex,:,:]
    plt.imshow(nf1)
    plt.colorbar()
    plt.show()

def get_int_box_data(field):
    # check if concatenated box exists
    f = match_file(cfg.data_dir,'{0}*lighttravel_cat.npy'.format(cfg.box_headers[field]))
    if f != None:  
        zi,zf,_,box_zMpc,box_xyMpc = [float(i) for i in re.findall("[-+]?\d+[\.]?\d*",f)] # regex magic that extracts numbers from a string
        print zi, zf, box_zMpc
        cfg.pms['zMpc'] = int(box_zMpc); cfg.pms['xyMpc'] = int(box_xyMpc)
        cfg.pms['zi'] = zi; cfg.pms['zf'] = zf
        box = np.load(cfg.data_dir+f)
        cfg.pms['shape']=box.shape
        return box
    # create and save concatenated box
    flist = match_files(cfg.data_dir,'{0}*lighttravel'.format(cfg.box_headers[field]))
    box_list = []; zi_list = []; zf_list = []; box_zMpc = 0;
    if len(flist)==0: raise RuntimeError('No file found for field {0}'.format(field))
    for f in flist:
        print f
        zi,zf,_,box_size,box_Mpc = [float(i) for i in re.findall("[-+]?\d+[\.]?\d*",f)] # regex magic that extracts numbers from a string
        print zi, zf, box_size, box_Mpc
        data_file = cfg.data_dir+f
        box_data = np.fromfile(data_file,dtype=np.float32)
        print "box_data.shape = ",box_data.shape
        #box_size = int(f.split('_')[-3])
        box_list.append(box_data.reshape((box_size,box_size,box_size)))
        zi_list.append(zi)
        zf_list.append(zf)
        box_zMpc += box_Mpc
    cfg.pms['zMpc'] = int(box_zMpc)
    cfg.pms['xyMpc'] = int(box_Mpc)
    cfg.pms['zi'] = min(zi_list)
    cfg.pms['zf'] = max(zf_list)
    sorted_box_list = [box for (zi,box) in sorted(zip(zi_list,box_list))]
    catbox = np.concatenate(sorted_box_list,axis=2)
    cfg.pms['shape']=catbox.shape
    print 'catbox',catbox.shape
    np.save('{5}{0}zstart{1}_zend{2}_FLIPBOXES1_{3}_{4}Mpc_lighttravel_cat'.format(cfg.box_headers[field],cfg.pms['zi'],cfg.pms['zf'],int(box_zMpc),int(box_xyMpc),cfg.data_dir),catbox)
    return catbox 

def plot_int_box_data(field):
    nf1 = get_int_box_data(field)
    chosenIndex=200
    print "nf1.shape = ",nf1.shape
    print "max = ",np.max(np.abs(nf1))
    print "min = ",np.min(np.abs(nf1))
    print "avg = ",np.mean(np.abs(nf1))
    nf1_slice = nf1[chosenIndex,:,:]
    plt.imshow(nf1_slice,cmap='jet')
    plt.colorbar()
    plt.show()

if __name__=='__main__':
  # nf = get_int_box_data('nf')
  # plot_box_data('nf',6.0)
  plot_int_box_data('vx')
