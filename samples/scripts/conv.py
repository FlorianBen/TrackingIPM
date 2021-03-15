#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

def gaussian(x, mu=0, sig=0.01):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


def main():
    x = np.linspace(-0.05, 0.05, 600)
    sig_pos = np.repeat([0., 1, 0.], 200)
    sig_gaus = gaussian(x)
    sing_conv = signal.convolve(sig_pos, sig_gaus, mode='same') / sum(sig_gaus)
    plt.plot(x, sig_pos)
    plt.plot(x, sig_gaus)
    plt.plot(x, sing_conv)
    plt.show()
    pass


if __name__ == "__main__":
    main()
