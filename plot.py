import numpy             as np
import matplotlib.pyplot as plt

bin_labels = np.loadtxt(fname="histogramdata", dtype=float, usecols=(0))
hist_data  = np.loadtxt(fname="histogramdata", dtype=float, usecols=(1))
bin_labels = [round(x,2) for x in bin_labels]

plt.plot(bin_labels, hist_data)
plt.show()
