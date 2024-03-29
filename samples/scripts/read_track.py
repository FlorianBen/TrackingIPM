import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def main():
    print('Main')
    f = h5py.File('build/tracking.h5', 'r')
    fig = plt.figure()

    bins = np.linspace(-0.025, 0.025, 100)

    countsi = []
    countsf = []

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for i in range(2000):
        pos_dset = np.asarray(f['track_' + str(i) + '/pos'][:])
        ax.plot(pos_dset['x'], pos_dset['y'], pos_dset['z'])
        countsi.append(pos_dset['y'][0])
        countsf.append(pos_dset['y'][-1])

    # plt.hist(countsi, bins, histtype='step')
    # plt.hist(countsf, bins, histtype='step')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_xlim3d(-0.06, 0.06)
    ax.set_ylim3d(-0.06, 0.06)
    ax.set_zlim3d(-0.06, 0.06)

    plt.show()


if __name__ == "__main__":
    main()
