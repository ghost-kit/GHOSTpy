from spacepy import pycdf as cdf
import numpy as np
import glob
import datetime
import matplotlib.pylab as plt

path = "./data/combined_30358/"

files = glob.glob(path+"*.cdf")

print ("Files: {}".format(files))

for f in files:
    cdf_data = cdf.CDF(f)
    epoch = cdf_data.get('Epoch')[:]
    magb = cdf_data.get('Magnitude')[:]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(epoch, magb)
    plt.show(fig)

