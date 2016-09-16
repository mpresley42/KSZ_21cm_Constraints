import unittest
from Sim21cm import *

DATA_DIR = '/Users/mpresley/Research/KSZ_21cm_Constraints/data/'
RUN_DIR = 'RandSeed_111_Sigma8_0.81577_h_0.68123_Omm_0.30404_Omb_0.04805_ns_0.96670_Rmfp_35.00_Tvir_60000.0_Squiggly_40.00_lnAs_3.06400'

class Test21cmSim(unittest.TestCase):
    def test_get_params(self):
        sim = Sim21cm(DATA_DIR,RUN_DIR)
        self.assertTrue(sim.pms['RandSeed']==111.)
        self.assertTrue(sim.pms['h']==0.68123)
        self.assertTrue(sim.pms['Rmfp']==35.)

    def test_set_constants(self):
        sim = Sim21cm(DATA_DIR,RUN_DIR)
        self.assertTrue(sim.pms['Tvir']==60000.)
        self.assertTrue(sim.pms['h']==0.68123)

    def test_setup_iboxes(self):
        

if __name__=='__main__':
    unittest.main()


