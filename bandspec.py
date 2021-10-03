import sys
import yaml
import math
import random
import numpy as np
import matplotlib.pyplot as plt

def normal_distribution(x, mean, sigma):
    return np.exp(-1*((x-mean)**2)/(2*(sigma**2)))/(math.sqrt(2*np.pi) * sigma)

def saveFile(bandspec):
    bandpsecFile = 'bandspec_'+str(int(bandset[0]))+'-'+str(int(bandset[1]))+'MHz'
    print "save to file: %s.npy" %(bandpsecFile)
    np.save(bandpsecFile, bandspec)

def bandspec(params):
    f1 = params['f1']
    f2 = params['f2']
    nchan = params['nchan']
    chan_freqs = np.linspace(f1, f2, num=nchan)
    band_center = (bandset[0] + bandset[1])/2.
    band = bandset[1] - bandset[0]
    
    # 2*(5sigma) = band, mean = band center
    mean1, sigma1 = band_center, band/2./3.
    bandspec = normal_distribution(chan_freqs, mean1, sigma1)
    bandspec = bandspec/max(bandspec)

    plot = False
    #plot = True 
    if (plot):
        plt.plot(chan_freqs, bandspec)
        plt.axvline(x=bandset[0], c='r')
        plt.axvline(x=bandset[1], c='r')
        plt.xlabel("Freq (MHz)")
        plt.ylabel("Weight")
        plt.grid()
        plt.savefig("bandspec.png")
        plt.show()

    return bandspec

 


if __name__=="__main__":

    # reading telescope setting from simulate file
    fastParam = "params/FAST_19Beam.params"
    with open(fastParam, 'r') as fp:
        params=yaml.load(fp.read(), Loader=yaml.FullLoader)
        
    global bandset
    bandset = [float(i) for i in sys.argv[1:]]
    
    print "Telescope info: %s MHz, %s MHz, %s" %(params['f1'], params['f2'], params['nchan'])
    print "Band Setting: ", bandset
    bandspec = bandspec(params)
    saveFile(bandspec)
