import numpy as np
import pylab as plt

# Haar wavelet's mother wavelet function
mHaar = lambda x: np.piecewise(x, [x<0, x>=0, x>=0.5, x>=1], [0,1,-1,0])
# General Haar Function
Haar = lambda n,k,x: 2.**(n/2.)*mHaar((2.**n)*x-k)
#Haar = lambda n,k,x: mHaar((2.**n)*x-k)

def plot_basic_wavelets():
    x = np.linspace(-5., 5., 100)
    for ii in range(3):
        plt.plot(x,Haar(ii,0,x),label=r"$\psi_{%d,0}(x)$"%(ii,))
    plt.ylim([-2,2])
    plt.xlabel(r"$x$",fontsize=18)
    plt.ylabel(r"$\psi_{n,k}(x)$",fontsize=18)
    plt.legend()
    plt.show()

    x = np.linspace(-5., 5., 100)
    for ii in range(3):
        plt.plot(x,Haar(0,ii,x),label=r"$\psi_{0,%d}(x)$"%(ii,))
    plt.ylim([-2,2])
    plt.xlabel(r"$x$",fontsize=18)
    plt.ylabel(r"$\psi_{n,k}(x)$",fontsize=18)
    plt.legend()
    plt.show()

def wavelet_decomp(x,f,nvals,kvals):
    Cnk = np.zeros([len(nvals),len(kvals)])
    for ii in range(len(nvals)):
        for jj in range(len(kvals)):
            nkHaar = Haar(nvals[ii],kvals[jj],x)
            Cnk[ii,jj] = np.trapz(f*nkHaar,x=x)
    return Cnk

def wavelet_recov(Cnk,nvals,kvals,x):
    f = np.zeros_like(x)
    for ii in range(len(nvals)):
        for jj in range(len(kvals)):
            nkHaar = Haar(nvals[ii],kvals[jj],x)
            f += Cnk[ii,jj]*nkHaar
    return f

def test_decomp():
    nmax = 10.; kmax = 10.
    dx = 0.1*2.**(-nmax)

    x = np.arange(-4.,6.,dx)
    f = (2./np.pi)*np.arctan(x-2.)

    nvals = np.arange(-nmax,nmax)
    kvals = np.arange(-kmax,kmax)
    Cnk = wavelet_decomp(x,f,nvals,kvals)
    plt.matshow(Cnk,origin='lower')
    plt.colorbar()
    plt.savefig('Cnk_matrix.pdf') 
    plt.show()

    plt.clf()
    fr = wavelet_recov(Cnk,nvals,kvals,x)
    plt.plot(x,f,label='Original')
    plt.plot(x,fr,label='Recovered')
    plt.legend(loc='lower right')
    plt.savefig('recovered_function.pdf')
    plt.show()

def build_intuition():
    nmax = 10.; kmax = 10.
    dx = 0.1*2.**(-nmax)
    x = np.arange(-10.,10.,dx)
    nvals = np.arange(-nmax,nmax)
    kvals = np.arange(-kmax,kmax)
    avals = np.arange(0.5,3.0,0.5)
    
    fig, ax = plt.subplots(len(avals),2, sharex='col')
    for ii,a in enumerate(avals):
        f = (2./np.pi)*np.arctan(a*x)
        Cnk = wavelet_decomp(x,f,nvals,kvals)
        fr = wavelet_recov(Cnk,nvals,kvals,x)
        ax[ii,0].plot(x,f,x,fr,label=r"$a = %d$"%a)
        ax[ii,0].set_ylabel(r"$f(x)$",fontsize=18)
        im = ax[ii,1].matshow(Cnk,origin='lower',vmin=-1,vmax=2)
        ax[ii,1].set_ylabel(r"$k$",fontsize=18)
        ax[ii,1].xaxis.set_ticks_position('bottom')
        ax[ii,1].set_aspect('equal',adjustable='box-forced')
        #ax[ii,0].set_aspect('equal',adjustable='box-forced')
    ax[-1,0].set_xlabel(r"$x$",fontsize=18)
    ax[-1,1].set_xlabel(r"$n$",fontsize=18)
    ax[0,0].set_title(r"$f(x)=\frac{2}{\pi}\arctan(ax)$",fontsize=18)
    ax[0,1].set_title(r"$C_{n,k}$",fontsize=18)

    cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    fig.colorbar(im, cax=cax)
    plt.savefig('build_intuition.pdf')
    #plt.show()
    

if __name__=='__main__':
    # plot_basic_wavelets()
    # test_decomp()
    build_intuition()
