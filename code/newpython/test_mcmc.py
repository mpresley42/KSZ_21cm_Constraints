import numpy as np
import pylab as plt
from mcmc import *

# Gaussian jump function
def jump_gauss(n=2,sig=0.1): 
    return np.random.normal(0.,sig,n)
# Linear model
def model_lin(x,p): 
    return p[0]*x+p[1]
# Guassian model
def model_gauss(x,p):
    return p[0]*np.exp(-(x-p[1])**2/p[2]**2)
# Generate data
def data_gen(x,model,truth,sigy):
    return model(x,truth)+np.random.normal(0.,sigy,len(x))

def test_plot_data_lin():
    # generate data
    p0 = [0.5,2.] # m,b
    lab = ['m','b']
    dx = np.arange(-5.,5.5,0.5)
    sigy = 0.5*np.ones_like(dx)
    dy = data_gen(dx,model_lin,p0,sigy)
    # instantiate mcmc class
    gm = GaussMCMC(dx,dy,np.diag(sigy),model_lin,jump_gauss,lab)
    # plot data
    gm.plot_data(truth=p0)

def test_plot_data_gauss():
    # generate data
    p0 = [5.,4.,3.] # A,mu,sig
    lab = ['A','mu','sig']
    dx = np.arange(-5.,15.,1.)
    sigy = 0.5*np.ones_like(dx)
    dy = data_gen(dx,model_gauss,p0,sigy)
    # instantiate mcmc class
    gm = GaussMCMC(dx,dy,np.diag(sigy),model_gauss,jump_gauss,lab)
    # plot data
    gm.plot_data(truth=p0)

def test_stair_lin(p0=[0.5,2.],p1=[1.,3.],nn=10000,lab=(r'$m$',r'$b$')):
    # generate data
    dx = np.arange(-5.,5.5,0.5)
    sigy = 0.5*np.ones_like(dx)
    dy = data_gen(dx,model_lin,p0,sigy)
    # instantiate mcmc class
    gm = GaussMCMC(dx,dy,np.diag(sigy),model_lin,jump_gauss,lab)
    # plot staircase
    gm.stair(p1,nn,truth=p0)

def test_stair_gauss(p0=[5.,4.,3.],p1=[3.,6.,2.],nn=10000,lab=(r'$A$',r'$\mu$',r'$\sigma$')):
    # generate data
    dx = np.arange(-5.,5.5,0.5)
    sigy = 0.5*np.ones_like(dx)
    dy = data_gen(dx,model_gauss,p0,sigy)
    # instantiate mcmc class
    gm = GaussMCMC(dx,dy,np.diag(sigy),model_gauss,jump_gauss,lab)
    # plot staircase
    gm.stair(p1,nn,truth=p0,fname='stair_test_blank.pdf')

if __name__=='__main__':
    #test_plot_data_gauss()
    #plt.show()
    test_stair_gauss(nn=1000)
    #test_stair_lin()
    plt.show()
