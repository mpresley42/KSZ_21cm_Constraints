import numpy as np
import pylab as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
import seaborn as sns

#mpl.rc('text', usetex = True)
sns.set_style("ticks",
    {'xtick.direction': u'in',
    'ytick.direction': u'in'})

class GaussMCMC:
    def __init__(self,dx,dy,cv,model,jump,lab):
        self.dx = dx # x data
        self.dy = dy # y data
        self.cv = cv # covariance matrix
        self.npm = len(lab) # number of parameters
        self.incv = self.__gen_inv()
        self.nd = len(dy) # number of data points
        self.model = model
        self.jump = jump
        # plot formatting
        self.lab = lab # labels for parameters
        self.cmap = sns.cubehelix_palette(light=1, as_cmap=True)
        self.clist = sns.color_palette("hls", self.npm)
        
    def __gen_inv(self):
        return np.linalg.pinv(self.cv)
        
    def like(self,p):
        diff = (self.dy-self.model(self.dx,p))
        return np.exp(-0.5*np.dot(np.transpose(diff),np.dot(self.incv,diff)))

    def __metro(self,p):
        delta_p = self.jump(self.npm)
        like_old = self.like(p)
        like_new = self.like(p+delta_p)
        if like_new>like_old:
            return delta_p
        elif np.random.uniform()<like_new/like_old:
            return delta_p
        else:
            return 0
    
    def chain(self,start,num):
        chain = np.zeros((num,len(start)))
        chain[0] = start
        for ii in xrange(1,num):
            chain[ii] = chain[ii-1]+self.__metro(chain[ii-1])
        return chain
 
    def stair(self,p0,num,truth=None,fname=None):
        # generate mcmc chain
        c1 = self.chain(p0,num)
        # clip initial points before convergence
        c1 = c1[int(0.1*num):,:]
        # set up figure
        #fig = plt.figure()
        #gs = gridspec.GridSpec(self.npm, self.npm)
        f, ax = plt.subplots(self.npm,self.npm,sharex=True,sharey=True)

        # loop over diagonal -- single variable histograms
        for ii in xrange(self.npm):
            #ax = fig.add_subplot(gs[ii,ii])
            #n, bins, patches = ax.hist(c1[:,ii],label=labels[ii],
            #    edgecolor=self.clist[ii],fill=False,histtype='step')
            sns.distplot(c1[:,ii],ax=ax[ii,ii],
                        label=self.lab[ii],color=self.clist[ii])
            if truth!=None: ax[ii,ii].axvline(x=truth[ii],c='k')
            ax[ii,ii].legend()
        # loop over stairs -- pair variables
        for ii in xrange(1,self.npm):
            for jj in xrange(0,ii):
                    #ax = fig.add_subplot(gs[ii,jj])
                    sns.kdeplot(c1[:,jj],c1[:,ii],ax=ax[ii,jj],
                            shade=True,shade_lowest=False,cmap=self.cmap)
                    if truth!=None:
                        ax[ii,jj].scatter((truth[jj],),(truth[ii],),c='k')
                        ax[ii,jj].axvline(x=truth[jj],c='k')
                        ax[ii,jj].axhline(y=truth[ii],c='k')
                    if jj==0: ax[ii,jj].set_ylabel(self.lab[ii])
                    if ii==self.npm-1:ax[ii,jj].set_xlabel(self.lab[jj])
        # turn off axes above the diagonal
        for ii in xrange(0,self.npm):
            for jj in xrange(ii+1,self.npm):
                ax[ii,jj].axis('off')
        # create plot of data and fit
        ax2 = f.add_axes([0.65, 0.65, 0.3, 0.3])
        ax2 = self.plot_data(truth=truth,ax=ax2,markersize=5)
        #sns.pointplot(x=self.dx, y=self.dy, ci=self.cv.diagonal(), 
        #            capsize=.2, ax=ax2)
        if truth!=None: 
            ax2.plot(self.dx,self.model(self.dx,truth),c='orange')
        if fname!=None: plt.savefig(fname)

    def plot_data(self,truth=None,ax=None,**kwargs):
        if ax==None: f,ax = plt.subplots(1,1)
        (_, caps, _) = ax.errorbar(self.dx, self.dy, yerr=self.cv.diagonal(), 
                    fmt='o',c='steelblue',**kwargs)
        for cap in caps:
            cap.set_markeredgewidth(1)
        if truth!=None:
            plt.plot(self.dx,self.model(self.dx,truth),c='orange')
        return ax




