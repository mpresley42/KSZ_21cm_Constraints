import sys
import os
import cfg
import numpy as np 
import load_data as ld
from redshift import *
import pylab as plt

def get_i_box_data(field,ibox,itot=8):
    # check if ibox exists
    f = match_file(cfg.data_dir,'{0}*lighttravel_cat_{1}_*.npy'.format(cfg.box_headers[field],ibox))
    if f != None:  
        zi,zf,_,box_zMpc,box_xyMpc,_,itot = [float(i) for i in re.findall("[-+]?\d+[\.]?\d*",f)] # regex magic that extracts numbers from a string
        print zi, zf, box_zMpc
        cfg.pms['zMpc'] = int(box_zMpc)/itot; cfg.pms['xyMpc'] = int(box_xyMpc)
        cfg.pms['itot'] = int(itot)
        cfg.pms['zi'] = zi; cfg.pms['zf'] = zf
        box = np.load(cfg.data_dir+f)
        cfg.pms['shape']=box.shape
        return box
    # create and save separated boxes 
    cat_box = ld.get_int_box_data(field)
    ibox_list = np.array_split(cat_box,itot,axis=2)
    cat_box=None
    for ii,box in enumerate(ibox_list):
        np.save('{5}{0}zstart{1}_zend{2}_FLIPBOXES1_{3}_{4}Mpc_lighttravel_cat_{6}_{7}'.format(
            cfg.box_headers[field],cfg.pms['zi'],cfg.pms['zf'],cfg.pms['zMpc'],cfg.pms['xyMpc'],cfg.data_dir,ii,itot),box)
    cfg.pms['itot'] = itot
    cfg.pms['zMpc'] = cfg.pms['zMpc']/itot
    cfg.pms['shape'] = ibox_list[ibox].shape
    return ibox_list[ibox]

def get_i_z_d(ibox,proper=False):
    z,d = get_z_d(cfg.pms['zi'],cfg.pms['zf'],
                    dlen=cfg.pms['shape'][2]*cfg.pms['itot'],proper=proper)
    zlist = np.array_split(z,cfg.pms['itot'])
    dlist = np.array_split(d,cfg.pms['itot'])
    return zlist[ibox],dlist[ibox]

def compute_i_tau(ibox,pretty=False):
    nf = get_i_box_data('nf',ibox)
    density = get_i_box_data('density',ibox)
    
    if ibox==0:
        z1,d1 = get_i_z_d(ibox,proper=True)
        z0,d0 = get_z_d(0,z1[0],proper=True)
        z = np.concatenate((z0,z1))
        d = np.concatenate((d0,d1))
    else:
        z,d = get_i_z_d(ibox,proper=True)

    nb = cfg.pms['nb0']*(1+z)**3
    chiHII = np.ones_like(z)
    if ibox==0: chiHII[len(d0):] = np.average((1.-nf)*(1+density),axis=(0,1))
    else:       chiHII = np.average((1.-nf)*(1+density),axis=(0,1))
    chiHeIII = np.zeros_like(chiHII)
    chiHeIII[np.where(z<3)] = 1

    fn = cfg.pms['sigT']*nb*(chiHII+0.25*chiHeIII*cfg.pms['YBBNp'])*mToMpc
    tau = sp.integrate.cumtrapz(fn,x=d,initial=0)
    
    # NOTE: This returns the proper distance!
    if ibox==0: return z1,d1,tau[len(d0):]
    else:       return z,d,tau

def compute_tot_tau(itot):
    tau0 = 0.
    tau_list = []; zlist = []; dplist = []
    for ii in range(itot):
        z,dp,tau = compute_i_tau(ii,pretty=False)
        tau_list.append(tau+tau0)
        zlist.append(z); dplist.append(dp)
        tau0 = tau[-1]+tau0
    tot_tau = np.concatenate(tau_list)
    tot_z = np.concatenate(zlist)
    tot_dp = np.concatenate(dplist)
    return tot_z,tot_dp,tot_tau

def compute_i_ndotq(ibox):
    nf = get_i_box_data('nf',ibox)
    density = get_i_box_data('density',ibox)
    # nf, density
 # Find the density and ionization weighting for the velocity field q
    q0 = (1-nf)*(1+density)
    nf = None; density = None
    #np.save('{0}/old_generated_data/q0.npy'.format(cfg.data_dir),q0)
    # q0
    z,d = get_i_z_d(ibox)
 # Get the spatial coordinates of all box cells
    xyzgd = get_i_xyz_gd(ibox)
    rgd = np.sqrt(xyzgd[0]*xyzgd[0]+xyzgd[1]*xyzgd[1]+xyzgd[2]*xyzgd[2])
    # q0, xyzgd (x3), rgd
 # Find the peculiar velocity 
    ndotq = np.zeros_like(q0)
    for ii,vii in enumerate(('vx','vy','vz')):
        vi = get_i_box_data(vii,ibox) #NOTE: These are co-moving velocities
        ndotq += q0*vi*xyzgd[ii]/rgd
    # q0, xyzgd (x3), rgd, ndotq, vi
    vi = None; q0=None; xyzgd=None; rgd=None
    if True: np.save('{0}/old_generated_data/ndotq_{1}'.format(cfg.data_dir,ibox),ndotq)
    return ndotq

# Create a grid of x,y,z coordinates for the standard box
def get_i_xyz_cd(ibox):
    box = cfg.pms['shape']
    boxMpc = np.array([cfg.pms['xyMpc'],cfg.pms['xyMpc'],cfg.pms['zMpc']])
    xcd = np.arange(-boxMpc[0]/2,boxMpc[0]/2,float(boxMpc[0])/box[0])
    ycd = np.arange(-boxMpc[1]/2,boxMpc[1]/2,float(boxMpc[1])/box[1])
    z,d = get_i_z_d(ibox)
    zcd = d[0] + np.arange(0,boxMpc[2],float(boxMpc[2])/box[2])
    return xcd,ycd,zcd
def get_i_xyz_gd(ibox):
    xcd,ycd,zcd = get_i_xyz_cd(ibox)
    xyzgd = np.meshgrid(xcd,ycd,zcd)
    return xyzgd

# Create 1D arrays of nx,ny,r coords for the standard box
def get_i_nxnyr_cd(ibox):
    box = cfg.pms['shape']
    boxMpc = np.array([cfg.pms['xyMpc'],cfg.pms['xyMpc'],cfg.pms['zMpc']])
    lx = boxMpc[0]/2.; ly = boxMpc[1]/2.; lz = boxMpc[2]
    z,d = get_i_z_d(ibox)
    
    # back of box -- throws away half the box but whatever
    df = d[0]+lz
    nx_max = lx / np.sqrt(lx*lx+df*df) # nx_min = - nx_max
    ny_max = ly / np.sqrt(ly*ly+df*df) # ny_min = - ny_max
    r_max = df
    r_min = np.sqrt(d[0]*d[0]+lx*lx+ly*ly) # r_min = d[0]

    nxcd = np.linspace(-nx_max,nx_max,box[0])
    nycd = np.linspace(-ny_max,ny_max,box[1])
    rcd = np.linspace(r_min,r_max,box[2])
    return nxcd,nycd,rcd

def compute_i_kSZ_linear(ibox,ndotq=None):
    if ndotq==None: ndotq = compute_i_ndotq(ibox)
    print "Have ndotq!"
    # ndotq    
 # Define box shape parameters
    box = cfg.pms['shape']
    boxMpc = np.array([cfg.pms['xyMpc'],cfg.pms['xyMpc'],cfg.pms['zMpc']])
    lx = boxMpc[0]/2.; ly = boxMpc[1]/2.; lz = boxMpc[2]*cfg.pms['itot']
    dxyz = boxMpc/box
 # Get tau
    zred,dp,tau = compute_i_tau(ibox)
    zred,d = get_i_z_d(ibox,proper=False)
    # ndotq, tau
    # Get n,r coord
    nxcd,nycd,rcd = get_i_nxnyr_cd(ibox)
    # Loop over angles
    dTkSZ = np.zeros([len(nxcd),len(nycd)])
    for ii in range(128):#len(nxcd)):
        print ii
        for jj in range(128):#len(nycd)):
            # Loop over z direction
            for kk in range(len(rcd)-1):
                # Find corresponding xyz coord
                x = nxcd[ii]*rcd[kk]
                y = nycd[jj]*rcd[kk]
                z = np.sqrt(rcd[kk]*rcd[kk]-x*x-y*y)
                # Find closest xyz indices
                xi = int(np.floor((lx+x)/dxyz[0]))
                yj = int(np.floor((ly+y)/dxyz[1]))
                zk = int(np.floor((z-d[0])/dxyz[2]))
                #dTkSZ[ii,jj] += np.exp(-tau[zk])*ndotq[xi,yj,zk]*(1+zred[zk])*(tau[zk+1]-tau[zk])
                dTkSZ[ii,jj] += np.exp(-tau[zk])*ndotq[xi,yj,zk]*(1+zred[zk])*(d[zk+1]-d[zk])
    dTkSZ = dTkSZ*cfg.pms['sigT']*cfg.pms['nb0']/cfg.pms['c']*mToMpc**2*cfg.pms['Tcmb']
    if True: np.save('{0}/old_generated_data/dTkSZ_{1}'.format(cfg.data_dir,ibox),dTkSZ)
    return dTkSZ 

if __name__=='__main__':
    ibox = int(sys.argv[1])
    nf = get_i_box_data('nf',ibox,itot=8)

    # xyzcd = get_i_xyz_cd(0)
    # np.savez('{0}/old_generated_data/xyzcd0.npz'.format(cfg.data_dir),*xyzcd)

    # for ii in range(8):
    #     z,d = get_i_z_d(ii)
    #     plt.plot(z,d)
    # plt.show()

    # tau0 = 0.
    # for ii in range(8):
    #     #z,d = get_i_z_d(ii)
    #     z,dp,tau = compute_i_tau(ii,pretty=False)
    #     plt.plot(z,tau+tau0)
    #     tau0 = tau[-1]+tau0
    # plt.xlabel(r'$z$',fontsize=18)
    # plt.ylabel(r'$\tau$',fontsize=18)
    # plt.show()
    # plt.savefig('{0}/figures/tau_all.pdf'.format(cfg.data_dir))
    
    # tot_z,tot_dp,tot_tau = compute_tot_tau(8)
    # plt.plot(tot_z,tot_tau)
    # plt.show()

    # if os.path.isfile('{0}ndotq_{1}.npy'.format(cfg.data_dir,ibox)):
    #     ndotq = np.load('{0}ndotq_{1}.npy'.format(cfg.data_dir,ibox))
    # else:
    #     ndotq = compute_i_ndotq(ibox) 

    #ndotq = compute_i_ndotq(ibox)

    ndotq = np.load('{0}/old_generated_data/ndotq_{1}.npy'.format(cfg.data_dir,ibox))
    dTkSZ = compute_i_kSZ_linear(ibox,ndotq)     
    
    # dTkSZ = compute_i_kSZ_test(ibox)
    # plt.imshow(dTkSZ,origin='lower')
    # plt.colorbar()
    # plt.savefig('{0}figures/vz_test_{1}.pdf'.format(cfg.data_dir,ibox))

    



