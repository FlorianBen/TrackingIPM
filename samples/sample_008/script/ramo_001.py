import numpy as np
import matplotlib.pyplot as plt
import h5py as h5


def main():
    f1 = h5.File('build/out.h5', mode='r')
    dset_pot = f1['pot']
    data_pot = np.array(dset_pot['x'])

    fig1, ax1 = plt.subplots(1, 1)

    ax1.imshow(data_pot[:, :, 15])
    plt.show()


if __name__ == "__main__":
    main()
