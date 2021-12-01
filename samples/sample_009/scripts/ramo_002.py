import h5py
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
from scipy.optimize import leastsq
from quantiphy import Quantity


weight_gaus = np.asarray((9, 5, 3, 2, 1.5, 1, 0.9,
                          0.8, 0.8, 0.8, 0.8, 0.9, 1, 1.5, 2, 3, 5, 9))

def main():
    def fitfunc(p, x): return p[0]*(1/(p[2]*np.sqrt(2*np.pi))
                                    )*np.exp(-0.5*((x-p[1])/p[2])**2)+p[3]

    def errfunc(p, x, y): return (y - fitfunc(p, x))

    f = h5py.File('build/tracking.h5', 'r')
    x_bins = np.linspace(-30, 30, 100)
    x_fits = np.linspace(-30, 30, 1000)

    countsi = []
    countsf = []

    for i in range(1000):
        pos_dset = np.asarray(f['track_' + str(i).zfill(5) + '/pos'][:])
        #ax.plot(pos_dset['x'], pos_dset['y'], pos_dset['z'])
        countsi.append(pos_dset['x'][0])
        countsf.append(pos_dset['x'][-1])

    countsf = np.array(countsf)*1000
    countsi = np.array(countsi)*1000

    hist, bin = np.histogram(countsf, bins=x_bins)
    hist_init, bin_init = np.histogram(countsi, bins=x_bins)

    init = [3000, 0.5, 0.5, 0.5]

    out = leastsq(errfunc, init, args=(bin[:-1], hist))
    c = out[0]
    out_init = leastsq(errfunc, init, args=(bin[:-1], hist_init))
    c_init = out_init[0]

    f1 = h5.File('build/out.h5', mode='r')
    dset_pot = f1['pot']
    data_pot = np.array(dset_pot['x'])

    f2 = h5.File('build/current.h5', mode='r')

    n_strip = 18

    fig1, ax1 = plt.subplots(2, 1)
    #fig2, ax2 = plt.subplots(n_strip, 1)

    current = np.zeros((n_strip, 100))
    profile = np.zeros((n_strip, 1))
    x_profile_gauss = np.array((-20.54, -13.4, -8.26, -6.66, -4.79, -3.42, -2.35, -
                                1.38, -0.46, 0.46, 1.38, 2.35, 3.42, 4.79, 6.66, 8.26, 13.4, 20.54))
    x_profil_lin = np.linspace(-14.72, 14.72, 32)

    for part in list(f2.keys()):
        part = f2[part]
        dset_pos = np.array(part['pos'])
        dset_traj = np.array(part['traj'])
        #ax1[0].plot(dset_pos['x'], dset_pos['y'], 'r.')
        for idx, curr in enumerate(part['current'].keys()):
            dset_curr = np.array(part['current'][curr])
            t = np.linspace(0, 8.343048851293794e-09, 100)
            current[idx, :] = current[idx, :] + dset_curr
            #ax2[idx].plot(t, dset_curr)

    for i in range(0, n_strip):
        t = np.linspace(0, 8.343048851293794e-09, 100)
        profile[i] = np.trapz(current[i, :], x=t)/weight_gaus[i]

    plot_c = ['#0C5DA5', '#00B945', '#FF9500',
              '#FF2C00', '#845B97', '#474747', '#9e9e9e']

    ax1[0].imshow(np.sum(data_pot[:, :, :], axis=2))
    ax1[1].plot(x_profile_gauss, np.max(fitfunc(c, x_fits)) * profile/np.max(profile), 'r:')

    ax1[1].plot(bin_init[:-1], hist_init, marker='+',
             color=plot_c[0], linestyle='None', label='Initial distribution')
    ax1[1].plot(x_fits, fitfunc(c_init, x_fits),
             label='$\mu=${0:.3}\n$\sigma=${1:.3}'.format(Quantity(c_init[1]/1000, 'm'), Quantity(c_init[2]/1000, 'm')), color=plot_c[0])
    ax1[1].plot(bin[:-1], hist, marker='+',
             color=plot_c[1], linestyle='None', label='Final distribution')
    ax1[1].plot(x_fits, fitfunc(c, x_fits),
                label='$\mu=${0:.3}\n$\sigma=${1:.3}'.format(Quantity(c[1]/1000, 'm'), Quantity(c[2]/1000, 'm')), color=plot_c[1])

    #ax1.plot(dset_traj['x'], dset_traj['y'], '.w')

    #fig, ax2 = plt.subplots(1, 1)
    #ax2.plot(t, dset_curr)

    plt.show()


def print_name(name):
    print(name)


def print_attrs(name, obj):
    print(name)
    for key, val in obj.attrs.items():
        print("    %s: %s" % (key, val))


if __name__ == "__main__":
    main()
