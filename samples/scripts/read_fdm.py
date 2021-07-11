import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def main():
    print('Main')
    data = np.loadtxt('build/out.txt')
    nx = 600
    ny = 400
    im = np.reshape(data, (ny, nx))
    #data = np.tile(im, 5)
    plt.imshow(im)
    plt.show()


if __name__ == "__main__":
    main()
