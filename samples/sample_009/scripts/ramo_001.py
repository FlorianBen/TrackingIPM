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

    n_strip = 32

    current = np.zeros((n_strip, 100))
    profile = np.zeros((32, 1))

    for part in list(f2.keys()):
        part = f2[part]
        for idx, curr in enumerate(part['current'].keys()):
            dset_curr = np.array(part['current'][curr])
            t = np.linspace(0, 8.343048851293794e-09, 100)
            current[idx, :] = current[idx, :] + dset_curr

    t = np.linspace(0, 8.343048851293794e-09, 100)
    for i in range(n_strip):
        profile[i] = np.trapz(current[i, :], x=t)

    fig1, ax1 = plt.subplots(2, 1)
    fig2, ax2 = plt.subplots(n_strip, 1, figsize=(10, 90))

    ax1[0].imshow(np.sum(data_pot[:, :, :], axis=2))
    ax1[1].plot(profile)

    for i in range(n_strip):
        ax2[i].plot(current[i, :])



    # ax1.plot(dset_traj['x'], dset_traj['y'], '.w')

    # fig, ax2 = plt.subplots(1, 1)
    # ax2.plot(t, dset_curr)

    fig2.savefig('test.png')
    plt.show()


if __name__ == "__main__":
    main()
