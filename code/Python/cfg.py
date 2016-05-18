import numpy as np

run_dir='RandSeed_111_Sigma8_0.81577_h_0.68123_Omm_0.30404_Omb_0.04805_ns_0.96670_Rmfp_35.00_Tvir_60000.0_Squiggly_40.00_lnAs_3.06400'
data_dir='/Users/mpresley/Research/KSZ_21cm_Constraints/data/'+run_dir+'/'
#data_dir='/Users/mpresley/Research/KSZ_21cm_Constraints/data/small/'+run_dir+'/'

box_headers = {'density':"updated_smoothed_deltax_",'vx':"updated_vx_",'vy':"updated_vy_",'vz':"updated_vz_",'v':"updated_v_",'nf':"xH_nohalos_"} 

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

pms = get_params(run_dir)
pms = get_constants(pms)