#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5


def main():
    f = h5.File('build/write_fieldmap.h5', mode='r')
    dset = f['sequence']
    data_z = np.array(dset['Fz'])
    data_x = np.array(dset['Fx'])
    plt.plot(data_z[0, 0, :])
    #plt.quiver(data_x[200, :, :],data_z[200, :, :],scale=1000)
    plt.show()
    pass


if __name__ == "__main__":
    main()
