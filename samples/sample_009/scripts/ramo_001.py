import h5py
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import math


def main():
    f1 = h5.File('build/out.h5', mode='r')
    dset_pot = f1['pot']
    data_pot = np.array(dset_pot['x'])

    f2 = h5.File('build/current.h5', mode='r')

    n_strip = 18

    fig1, ax1 = plt.subplots(2, 1)
    fig2, ax2 = plt.subplots(n_strip, 1)

    current = np.zeros((n_strip, 100))
    profile = np.zeros((n_strip, 1))
    x_profile_gauss = np.array((-20.54, -13.4, -8.26, -6.66, -4.79, -3.42, -2.35, -
                         1.38, -0.46, 0.46, 1.38, 2.35, 3.42, 4.79, 6.66, 8.26, 13.4, 20.54))
    x_profil_lin = np.linspace(-14.72e-3, 14.72e-3, 32)

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
        profile[i] = np.trapz(current[i, :], x=t)

    ax1[0].imshow(np.sum(data_pot[:, :, :], axis=2))
    ax1[1].plot(x_profile_gauss,profile)

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
