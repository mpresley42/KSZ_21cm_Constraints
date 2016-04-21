import os
import subprocess
from collections import OrderedDict

# RandSeed_111_Sigma8_0.83_h_0.67_Omm_0.32_Omb_0.022_ns_0.96_Rmfp_30.00_Tvir_10000.0_Squiggly_31.5_lnAs_3.06400   
params = OrderedDict([('RandSeed',111),('Sigma8',0.83),('h',0.67),('Omm',0.32),
            ('Omb',0.022),('ns',0.96),('Rmfp',30.0),('Tvir',10000.0),
            ('Squiggly',31.5),('lnAs',3.064)])

#params = {'RandSeed':111,'Sigma8':0.83,'h':0.67,'Omm':0.32,'Omb':0.022,
#        'ns':0.96,'Rmfp':30.0,'Tvir':10000.0,'Squiggly':31.5,'lnAs':3.064}
base_dir = '/scratch1/scratchdirs/mpresley/21cm_FAST_Sims/grid_runs/'

def create_submit_script(label):
    # get text of template file
    with open('varyOne_FullBoxIonHist_template.sh', 'r') as file:
        file_text = file.readlines()

    # now change the lines for this run
    # workLoc="/scratch1/scratchdirs/mpresley/21cm_FAST_Sims/test6"
    file_text[29] = 'workLoc="{0}{1}"\n'.format(base_dir,label)

    # write everything to a new file
    run_script = 'varyOne_FullBoxIonHist_{0}.sh'.format(label)
    with open(run_script, 'w') as file:
        file.writelines(file_text)

    # submit the new script
    subprocess.call(["sbatch {0}".format(run_script),], shell=True)

def create_run_directory(**kwargs):
    # replace any default parameters that we are changing
    for var in kwargs:
        params[var] = kwargs[var]   
    # create a label for the run directory
    label = ''
    for var in kwargs:
        if label != '': label += '_'
        label += str(var) + '_' + str(kwargs[var])
    # make run directory
    run_dir = base_dir + label
    if not os.path.exists(run_dir):
        os.makedirs(run_dir)
    # create the varyTvir.dat file for the run directory
    file_text = ''
    for item in list(params.items()):
        if file_text != '': file_text += ' '
        #file_text += str(var) + ' ' + str(params[var])
        file_text += str(item[1])
    with open('{0}/varyTvir.dat'.format(run_dir), 'w') as file:
        file.write(file_text+'\n')
    return label

def do_run(**kwargs):
    label = create_run_directory(kwargs)
    create_submit_script(label)

if __name__=='__main__':
    #label = create_run_directory(Rmfp=20.)
    label = create_run_directory(Tvir=21000.)
    create_submit_script(label)

#     for Tvir in (10000., 12000., 14000.,):
#         do_run(Tvir=Tvir)
#     for Rmfp in (20., 30., 40.,):
#         do_run(Rmfp=Rmfp)
#     for Squiggly in (25., 30., 35.):
#         do_run(Squiggly=Squiggly)
