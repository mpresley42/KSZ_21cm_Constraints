import numpy as np 
import timeit
from Simbox import *

data_dir = '/Users/mpresley/Research/KSZ_21cm_Constraints/data/mesinger_1/original_1_mesinger_1/'

def test_kslice():
    start = timeit.default_timer()
    sim = Sim(data_dir)
    sim.box['density'].slice(100) 
    sim.box['density'].close()
    sim.box['density'].slice(100)
    end = timeit.default_timer()
    print "Time for one k-slice is: {0}".format(end-start)
    return end-start

def test_islice():
    start = timeit.default_timer()
    sim = Sim(data_dir) 
    sim.box['nf'].slicex(100)
    sim.box['density'].slicex(100)
    end = timeit.default_timer()
    print "Time for one i-slice is: {0}".format(end-start)
    return end-start

def test_jslice():
    start = timeit.default_timer()
    sim = Sim(data_dir) 
    sim.box['density'].slicey(100)
    end = timeit.default_timer()
    print "Time for one j-slice is: {0}".format(end-start)
    return end-start

def test_cube_read():
    start = timeit.default_timer()
    cube1 = np.memmap('/Users/mpresley/Research/KSZ_21cm_Constraints/data/mesinger_1/original_1_mesinger_1/updated_smoothed_deltax__zstart005.00000_zend009.56801_FLIPBOXES0_1024_1600Mpc_lighttravel')
    #cube1 = np.memmap('/Users/mpresley/Research/KSZ_21cm_Constraints/data/mesinger_1/original_1_mesinger_1/updated_smoothed_deltax__zstart005.00000_zend009.56801_FLIPBOXES0_1024_1600Mpc_lighttravel')
    cube[:,:,100]
    cube.close()
    end = timeit.default_timer()
    print "Time for one cube reading is: {0}".format(end-start)
    return end-start



if __name__=='__main__':
    # sim = Sim(data_dir)
    # for kk in xrange(1024):
    #     start = timeit.default_timer()
    #     # chiHII[kk+len(d0)] = np.average((1.-sim.box['nf'].slice(kk))* \
    #     #     (1.+sim.box['density'].slice(kk)))
    #     nf = sim.box['nf'].slice(kk)
    #     density = sim.box['density'].slice(kk)
    #     end = timeit.default_timer()
    #     print kk, end-start
    test_islice()