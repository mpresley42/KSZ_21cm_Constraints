import sys
import os
import subprocess
import fnmatch as fnm
import re

def match_files(location,pattern):
    files = []
    for f in os.listdir(location):
        if fnm.fnmatch(f, pattern):
            files.append(f)
    return files

# xH_nohalos_z005.50_nf0.000000_eff30.0_HIIfilter1_Mmin2.7e+09_RHIImax35_400_200Mpc
# delta_T_v3_no_halos_z005.50_nf0.000000_useTs0_zetaX-1.0e+00_alphaX-1.0_TvirminX-1.0e+00_aveTb000.00_Pop-1_400_200Mpc
def scrape_data(location):
    pattern = 'delta_T_v3_*Mpc'
    filenames = match_files(location,pattern)
    wf = open('{0}z_nf_Tb.dat'.format(location),'w')
    for f in filenames:
        nums = [float(i) for i in re.findall("[-+]?\d+[\.]?\d*",f)] # regex magic that extracts numbers from a string
        z = nums[1]
        nf = nums[2]
        Tb = nums[9]
        wf.write('{0} {1} {2}\n'.format(z,nf,Tb))
    wf.close()

if __name__=='__main__':
    #loc = '/scratch1/scratchdirs/mpresley/21cm_FAST_Sims/test6/RandSeed_111_Sigma8_0.81577_h_0.68123_Omm_0.30404_Omb_0.04805_ns_0.96670_Rmfp_35.00_Tvir_60000.0_Squiggly_30.00_lnAs_3.06400/Boxes/'
    #loc = './'
    if len(sys.argv)==1:
        loc = './'
    else:
        loc = sys.argv[1]
    scrape_data(loc)
