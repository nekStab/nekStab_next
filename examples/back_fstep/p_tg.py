#!/usr/bin/env python
import matplotlib.pyplot as plt
from numpy import logspace
import numpy as np
params = {'text.usetex': False,
          'font.size': 8,
          'legend.fontsize': 8,
          'legend.handlelength': 2.5,}
plt.rcParams.update(params)
plt.style.use('seaborn-white')

formt = 'png'
ajust = 'tight'
qual = 500
fig_width = 3.5
fig_height = 2.45

class res(object):
    def __init__(self, filename):
        #data = np.transpose(np.loadtxt(filename))
        print('Reading '+filename)
        data = np.transpose(np.genfromtxt(filename))
        self.t  = data[0]
        self.e  = data[1]
        del data

class SpectreS(object):
    def __init__(self, filename):
        print('Reading '+filename)
        data = np.transpose(np.genfromtxt(filename))
        self.sigma  = data[0]
        self.res   = data[1]
        del data

########################################

if __name__ == '__main__':

    fig=plt.figure()
    fig.set_size_inches(fig_width, fig_height)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$E(t)/E(0)$')

    f = res('barkley2008_fig5.ref')
    plt.plot(f.t,f.e,c='k',lw=0.2,ls='-',label='Barkley et al. (2008)')

    base = "transient_growth"  # reference folder to copy
    p1 = logspace(0, 2, 10)
    for i in range(len(p1)):
        try:
            p1[i] = round(p1[i],2)
            path = "t_" + str(p1[i])
            f = SpectreS(path+'/Spectre_Sp.dat')
            print(p1[i],f.sigma[0])
            plt.scatter(p1[i],f.sigma[0]**2,s=3, facecolors='none', edgecolors='k')
        except:
            print(' File not found...')
            pass
    plt.legend(loc='best',fontsize=6)
    fname='transient_growth_Re500.'+formt
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
    print('------------------------------------------')
