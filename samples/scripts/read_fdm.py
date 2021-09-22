import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from matplotlib.colors import LogNorm


def dist(pos1, pos2):
    pos_t = pos1 - pos2
    dist = np.sqrt(pos_t[0]**2+pos_t[1]**2)
    return dist


def main():
    print('Main')
    data = np.loadtxt('build/out.txt')
    nx = 2000
    ny = 250
    dx = 0.08/2000

    im = np.reshape(data, (ny, nx))
    # E = np.gradient(im, dx)

    # v0 = np.array([1e3, 1.e4])
    # time = np.linspace(0, 1000e-9, 100)

    # pos = np.array([v0*t for t in time])

    # x, y = np.linspace(-0.04, 0.04, nx), np.linspace(0, 0.01, ny)
    # xi, yi = np.meshgrid(x, y)
    # XY = np.dstack((xi, yi))

    # for j in range(np.size(y)):
    #     for i in range(np.size(x)):
    #         dist(pos, XY[j,i])

    fig, axs = plt.subplots(1, 1)
    # axs[0].imshow(E[0])
    # axs[1].imshow(E[1])
    # axs.imshow(10*np.log(np.sqrt(E[0]**2+E[1]**2)),
    #            extent=(-0.04, 0.04, 0, 0.01))
    # axs.plot(pos[:, 0], pos[:, 1], 'r')

    axs.imshow(im)

    plt.show()


if __name__ == "__main__":
    main()
