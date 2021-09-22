import h5py
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import math


def main():
    f = h5.File('build/ramo.h5', mode='r')
    dset = f['pot']
    data_x = np.array(dset['x'])

    dset_pos = np.array(f['pos'])
    dset_traj = np.array(f['traj'])
    dset_curr = np.array(f['current/current_9'])

    t =  np.linspace(0, 8.343048851293794e-09, 100)
    print(np.trapz(dset_curr, x=t))

    fig, ax1 = plt.subplots(1, 1)
    ax1.imshow(data_x[:, :, 9])
    ax1.plot(dset_pos['x'], dset_pos['y'], 'rx')
    #axs.plot(dset_traj['x'], dset_traj['y'], '.w')

    fig, ax2 = plt.subplots(1, 1)
    ax2.plot(t, dset_curr)

    plt.show()


if __name__ == "__main__":
    main()
