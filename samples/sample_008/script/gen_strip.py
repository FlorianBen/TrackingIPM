#!/usr/bin/env python
import numpy as np


def main():
    """
    Fonction principale.
    :return:
    """

    inter = 0.32
    strip_size = 2.5
    n_strips = 32

    bins = np.zeros(n_strips * 2)
    bins[::2] += inter
    bins[1::2] += strip_size
    bins = np.insert(np.cumsum(bins), 0, 0)
    bins = bins - np.max(bins) / 2 - inter/2
    bins = bins[1:]
    bins = np.reshape(bins, (n_strips, 2))
    np.savetxt('my_strips.txt', bins, fmt='%.2f')


if __name__ == '__main__':
    main()
