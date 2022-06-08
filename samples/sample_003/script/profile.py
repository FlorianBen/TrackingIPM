#!/usr/bin/env python
from operator import pos
import h5py
import numpy as np
import matplotlib.pyplot as plt


def main():
    """
    Fonction principale.
    :return:
    """
    f = h5py.File('build/tracking.h5', 'r')

    fig = plt.figure()
    ax = fig.add_subplot(111)

    posi_dset = np.asarray(f['track/posi'][:])
    posf_dset = np.asarray(f['track/posf'][:])

    x_bins = np.linspace(-25, 25, 50)*1e-3
    x_fits = np.linspace(-25, 25, 1000)*1e-3

    hist, bin = np.histogram(posi_dset['x'], bins=x_bins)
    hist_f, bin = np.histogram(posf_dset['x'], bins=x_bins)

    ax.plot(bin[:-1], hist, marker='+', linestyle='None',
            label='Initial distribution')
    ax.plot(bin[:-1], hist_f, marker='x', linestyle='None',
            label='Final distribution')
    
    ax.legend()
    #ax.plot(posi_dset['x'])

    plt.show()


if __name__ == '__main__':
    main()
