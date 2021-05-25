#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5


def main():
    f = h5.File('build/write_fieldmap.h5', mode='r')
    dset = f['sequence']
    data_z = np.array(dset['z'])
    data_x = np.array(dset['x'])
    data_y = np.array(dset['y'])

    #plt.plot(data_x[0, :, 0])
    plt.quiver(data_x[::20, ::20, 0], data_y[::20, ::20, 0], scale=5e5)
    plt.show()
    pass


if __name__ == "__main__":
    main()
