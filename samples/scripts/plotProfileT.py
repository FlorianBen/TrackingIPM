import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
from scipy.optimize import leastsq
from quantiphy import Quantity
import h5py


def main():
    def fitfunc(p, x): return p[0]*(1/(p[2]*np.sqrt(2*np.pi))
                                    )*np.exp(-0.5*((x-p[1])/p[2])**2)+p[3]

    def errfunc(p, x, y): return (y - fitfunc(p, x))

     
    f = h5py.File('build/tracking.h5', 'r')
    x_bins = np.linspace(-25, 25, 100)
    x_fits = np.linspace(-25, 25, 1000)


    countsi = []
    countsf = []

    for i in range(20000):
        pos_dset = np.asarray(f['track_' + str(i) + '/pos'][:])
        #ax.plot(pos_dset['x'], pos_dset['y'], pos_dset['z'])
        countsi.append(pos_dset['y'][0])
        countsf.append(pos_dset['y'][-1])
        
        
    countsf = np.array(countsf)*1000
    countsi = np.array(countsi)*1000

    hist, bin = np.histogram(countsf, bins=x_bins)
    hist_init, bin_init = np.histogram(countsi, bins=x_bins)

    init = [3000, 0.5, 0.5, 0.5]

    out = leastsq(errfunc, init, args=(bin[:-1], hist))
    c = out[0]
    out_init = leastsq(errfunc, init, args=(bin[:-1], hist_init))
    c_init = out_init[0]

    plot_c = ['#0C5DA5', '#00B945', '#FF9500',
              '#FF2C00', '#845B97', '#474747', '#9e9e9e']

    scale = 1.5
    fig, ax1 = plt.subplots(figsize=[scale*5.40, scale*3.55])
    ax1.set_title('Profile measurement with particle tracking')
    ax1.set_xlabel('Transversal position ($\mathrm{mm}$)')
    ax1.set_ylabel('Counts ($\mathrm{A.U.}$)')
    ax1.plot(bin_init[:-1], hist_init, marker='+',
             color=plot_c[0], linestyle='None', label='Initial distribution')
    ax1.plot(x_fits, fitfunc(c_init, x_fits),
             label='$\mu=${0:.3}\n$\sigma=${1:.3}'.format(Quantity(c_init[1]/1000, 'm'), Quantity(c_init[2]/1000, 'm')), color=plot_c[0])
    ax1.plot(bin[:-1], hist, marker='+',
             color=plot_c[1], linestyle='None', label='Final distribution')
    ax1.plot(x_fits, fitfunc(c, x_fits),
             label='$\mu=${0:.3}\n$\sigma=${1:.3}'.format(Quantity(c[1]/1000, 'm'), Quantity(c[2]/1000, 'm')), color=plot_c[1])
    ax1.legend(loc='upper left')
    plt.show()

if __name__ == '__main__':
    plt.style.use('science')
    main()
