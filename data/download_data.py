import numpy as np
import subprocess as sp
import os
devnull = open(os.devnull, 'w')

# define remote and local data directories
sim_dir='mpresley@edison.nersc.gov:/scratch1/scratchdirs/mpresley/21cm_FAST_Sims/test4/'
local_data_dir='/Users/mpresley/Research/KSZ_21cm_Constraints/data/'

box_headers = {'density':"updated_smoothed_deltax_",'vx':"updated_vx_",'vy':"updated_vy_",'vz':"updated_vz_",'nf':"xH_nohalos_"}
box_ends = {'density':"_400_200Mpc",'vx':"_400_200Mpc",'vy':"_400_200Mpc",'vz':"_400_200Mpc",'nf':"_nf*_400_200Mpc"} #0.000000_eff30.0_HIIfilter1_Mmin8.4e+09_RHIImax35

redshifts=range(6,9,1)
fields=box_headers.keys()
# redshifts=(6.,)
# fields=('nf',)

# create a temporary directory to store the parameter file
temp_dir=local_data_dir+'temp/'
if not os.path.exists(temp_dir):
    sp.call(['mkdir',temp_dir])
# copy over the parameter file
param_file='varyTvir.dat'
sp.call(['scp','-r',sim_dir+param_file,temp_dir],stdout=devnull, stderr=devnull)
# create the directory name from the parameters
params = np.loadtxt(temp_dir+param_file)
# 'RandSeed_111_Sigma8_0.81577_h_0.68123_Omm_0.30404_Omb_0.04805_ns_0.96670_Rmfp_35.00_Tvir_60000.0_Squiggly_30.00_lnAs_3.06400'
run_dir='RandSeed_{0:g}_Sigma8_{1:01.5f}_h_{2:01.5f}_Omm_{3:01.5f}_Omb_{4:01.5f}_ns_{5:01.5f}_Rmfp_{6:02.2f}_Tvir_{7:05.1f}_Squiggly_{8:02.2f}_lnAs_{9:01.5f}'.format(*params)

data_dir='{0}{1}/Boxes'.format(sim_dir,run_dir)
local_run_dir='/Users/mpresley/Research/KSZ_21cm_Constraints/data/'+run_dir
if not os.path.exists(local_run_dir):
    sp.call(['mkdir',local_run_dir])

# # download field boxes
# for z in redshifts:
#     for field in fields:
#         data_file = '{0}/{1}z{2:06.2f}{3}'.format(data_dir,box_headers[field],z,box_ends[field])
#         sp.call(['scp','-r',data_file,local_run_dir],stdout=devnull, stderr=devnull)
#         print z, " ", field, " downloaded"
# sp.call(['cp',temp_dir+param_file,local_run_dir],stdout=devnull, stderr=devnull)

# # download mean field data files
# for field in ("avTb","nf","weightedxe"):
#     data_file = '{0}/{1}_*'.format(data_dir,field)
#     sp.call(['scp','-r',data_file,local_run_dir],stdout=devnull, stderr=devnull)

# download interpolated boxes
print data_dir+'/*lighttravel'
sp.call(['scp','-r',data_dir+'/*lighttravel',local_run_dir],stdout=devnull, stderr=devnull)
