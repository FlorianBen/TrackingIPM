import h5py
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import math


def dist(pos1, pos2):
    pos_t = pos1 - pos2
    dist = np.sqrt(pos_t[0]**2+pos_t[1]**2)
    return dist


def main():
    f = h5.File('build/ramo.h5', mode='r')
    dset = f['pot']
    data_x = np.array(dset['x'])

    dset_pos = np.array(f['pos'])
    print(dset_pos['x'])

    fig, axs = plt.subplots(1, 1)

    print(np.max(data_x[:, :, 0]))

    axs.imshow(data_x[:, :, 0])
    axs.plot(dset_pos['y'], dset_pos['x'], 'rx')

    plt.show()


if __name__ == "__main__":
    main()
