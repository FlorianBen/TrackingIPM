#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5


def main():
    f = h5.File('build/write_fieldmap.h5', mode='r')
    for i in range(8):
        dset = f['Efield{}'.format(i)]
        data_z = np.array(dset['z'])
        data_x = np.array(dset['x'])
        data_y = np.array(dset['y'])

        n = 100

        fig, axs = plt.subplots(3, 3)
        axs[0,0].imshow(data_x[:, :, n], clim=(-6e4, 6e4))
        axs[1,0].imshow(data_x[:, n, :])
        axs[2,0].imshow(data_x[n, :, :], clim=(-6e4, 6e4))
        axs[0,1].imshow(data_y[:, :, n])
        axs[1,1].imshow(data_y[:, n, :])
        axs[2,1].imshow(data_y[n, :, :])
        axs[0,2].imshow(data_z[:, :, n], clim=(-6e4, 6e4))
        axs[1,2].imshow(data_z[:, n, :])
        axs[2,2].imshow(data_z[n, :, :])


        axs[0,0].set_title('$E_{x}$ XY plane')
        axs[1,0].set_title('$E_{x}$ XZ plane')
        axs[2,0].set_title('$E_{x}$ YZ plane')
        axs[0,1].set_title('$E_{y}$ XY plane')
        axs[1,1].set_title('$E_{y}$ XZ plane')
        axs[2,1].set_title('$E_{y}$ YZ plane')
        axs[0,2].set_title('$E_{z}$ XY plane')
        axs[1,2].set_title('$E_{z}$ XZ plane')
        axs[2,2].set_title('$E_{z}$ YZ plane')
        #plt.show()
        plt.savefig('test{}.png'.format(i))
    #plt.quiver(data_x[::20, ::20, n], data_y[::20, ::20, n], scale=5e5)
    #plt.show()
    pass


if __name__ == "__main__":
    main()
