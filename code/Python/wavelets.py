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
    plt.imshow(Cnk,origin='lower')
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

if __name__=='__main__':
    plot_basic_wavelets()
    test_decomp()

