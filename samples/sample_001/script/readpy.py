import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib
matplotlib.use('TkAgg')


def main():
  print('Main')
  f = h5py.File('build/samples/sample_001/write_field.h5', 'r')
  dset = f['Field']
  data = np.reshape(dset, (-1, 400))

  time = f['Scalar/time']

  plt.plot(data[0]['x'])

  fig, ax = plt.subplots()
  xdata, ydata = [], []
  ln, = plt.plot([], [], '-')

  def init():
    ax.set_ylim(-1.1e5, 1.1e5)
    ax.set_xlim(-50, 50)
    return ln,

  def update(frame):
    ax.set_title('t: ' + str(time[frame]*1e9))
    xdata = np.linspace(-50, 50, 400)
    ydata = data[frame]['x']
    ln.set_data(xdata, ydata)
    return ln,

  ani = FuncAnimation(fig, update, frames=range(0, 1600),
                      init_func=init, blit=False, interval=30)
  ani.save("movie.mp4")


if __name__ == "__main__":
  main()
