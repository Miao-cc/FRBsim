import sys
import yaml
import math
import random
import numpy as np
import matplotlib.pyplot as plt


def normal_distribution(x, mean, sigma):
    return np.exp(-1*((x-mean)**2)/(2*(sigma**2)))/(math.sqrt(2*np.pi) * sigma)

# reading telescope setting from simulate file
fastParam = "simBinarydata/FAST.params"
with open(fastParam, 'r') as fp:
    params=yaml.load(fp.read())

bandset = [float(i) for i in sys.argv[1:]]

f1 = params['f1']
f2 = params['f2']
nchan = params['nchan']
print "Telescope info: %s MHz, %s MHz, %s" %(f1, f2, nchan)
print "Band Setting: ", bandset

chan_freqs = np.linspace(f1, f2, num=nchan)
#print chan_freqs[0:3]
band_center = (bandset[0] + bandset[1])/2.
band = bandset[1] - bandset[0]

# 2*(5sigma) = band, mean = band center
mean1, sigma1 = band_center, band/2./3.
bandspec = normal_distribution(chan_freqs, mean1, sigma1)
bandspec = bandspec/max(bandspec)

# if True, plt the distribution
#plot = True
plot = False
if (plot):
    plt.plot(chan_freqs, bandspec)
    plt.axvline(x=bandset[0], c='r')
    plt.axvline(x=bandset[1], c='r')
    plt.xlabel("Freq (MHz)")
    plt.ylabel("Weight")
    plt.grid()
    plt.show()
 
bandpsecFile = 'bandspec_'+str(int(bandset[0]))+'-'+str(int(bandset[1]))+'MHz'
print "save to file: ", bandpsecFile, ".npy"
np.save(bandpsecFile, bandspec)
