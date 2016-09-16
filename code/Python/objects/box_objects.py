
BOX_HEADERS = {'density':"updated_smoothed_deltax_",'vx':"updated_vx_",'vy':"updated_vy_",'vz':"updated_vz_",'v':"updated_v_",'nf':"xH_nohalos_"} 

class Params(dict):
    def __init__(self,run_dir):
        dict.__init__(self)
        self.set_constants(self)
        self.get_params(self,run_dir)

    # return a dictionary of parameters based off of the directory name
    def get_params(self,run_dir):
        pms=run_dir.split('_')
        for ii in xrange(0,len(pms),2):
            self[pms[ii]] = float(pms[ii+1])

    def set_constants(self):
        self['Sigma8']=0.83
        self['h']=0.67
        self['Omm']=0.32
        self['Omb']=0.022/(pms['h']*pms['h'])
        self['YBBNp']=0.24
        # =============================
        self['H0']=self['h']*3.241e-18 # s^-1 # = h * 100 km/s/Mpc
        self['mp']=1.672622e-27 # kg  
        self['G']=6.67384e-11 # N m^2 kg^-2
        self['c']=2.99792458e8 # m/s
        self['mHe']=6.6464764e-27 # kg
        self['mH']=1.6737236e-27 # kg
        self['mp']=1.6726219e-27 # kg
        self['mu']=1.0 + self['YBBNp']*0.25*(self['mHe']/self['mH']-1.0)
        self['nb0']=(3.*self['H0']**2*self['Omb'])/(8.*np.pi*self['G']*self['mu']*self['mp']) # at z=0. scales with (1+z)^3
        self['sigT']=6.65246e-29 # m^2
        self['Tcmb']=2.725e6 # microK

    def get_box_params(self,f):
        zi,zf,zMpc,xyMpc,itot = parse_ibox_filename(f)
        self.pms['zMpc'] = zMpc; self.pms['xyMpc'] = xyMpc
        self.pms['itot'] = itot
        self.pms['zi'] = zi; self.pms['zf'] = zf

f = match_file(self.data_dir,'{0}*lighttravel_cat_{1}_*.npy'.format(BOX_HEADERS[field],ibox))

def parse_ibox_filename(f):
    zi,zf,_,box_zMpc,box_xyMpc,_,itot = [float(i) for i in re.findall("[-+]?\d+[\.]?\d*",f)] # regex magic that extracts numbers from a string
    itot = int(itot)
    zMpc = int(box_zMpc)/itot; xyMpc = int(box_xyMpc)
    return zi,zf,zMpc,xyMpc,itot
    
class Box:
    def __init__(self,field,ibox,itot=8,params):
        self.pms = params
        self.box = get_i_box_data(field,ibox,itot)

    def get_i_box_data(self,field,ibox,itot=8):
        # check if ibox exists
        f = match_file(self.data_dir,'{0}*lighttravel_cat_{1}_*.npy'.format(BOX_HEADERS[field],ibox))
        if f != None:  
            zi,zf,zMpc,xyMpc,itot = parse_file(f)
            assert self.pms['zi'] == zi
            assert self.pms['zf'] == zf
            assert self.pms['zMpc'] == zMpc
            assert self.pms['xyMpc'] == xyMpc
            assert self.pms['itot'] == itot
            box = np.load(self.data_dir+f)
            return box
        # create and save separated boxes 
        cat_box = ld.get_int_box_data(field)
        ibox_list = np.array_split(cat_box,itot,axis=2)
        cat_box=None
        for ii,box in enumerate(ibox_list):
            np.save('{5}{0}zstart{1}_zend{2}_FLIPBOXES1_{3}_{4}Mpc_lighttravel_cat_{6}_{7}'.format(
                self.box_headers[field],self.pms['zi'],self.pms['zf'],self.pms['zMpc'],self.pms['xyMpc'],self.data_dir,ii,itot),box)
        # change to asserts
        assert self.pms['itot'] == itot
        assert self.pms['zMpc'] == self.pms['zMpc']/itot
        assert self.pms['shape'] == ibox_list[ibox].shape
        return ibox_list[ibox]






