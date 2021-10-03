import sys
import pywt
import yaml
import numpy as np

def usage(SettingFile):
    if len(SettingFile) > 0:
        print "Reading from file: ", SettingFile
    else: 
        print """usages:
        python simMultiPsr.py simulate.yaml
        -----------------------------------"""
        SettingFile = "simulate.yaml"
    return SettingFile

def inputCheck(fileList, fileLen):
    if len(fileList) != fileLen:
        print "Input Error: ", fileList, fileLen
        sys.exit(0)

def readBinary(filename, nchan, nsamp):
    """
    read simulation data
    need input filename, channel of the frequence(nchan), samples in time (nsmap)
    return numpy array in shape (nsamp, nchan)
    """
    with open(filename, 'rb') as binfile:
        # def read block
        # read block
        blockLen = nsamp*nchan
        binData = binfile.read(4*blockLen)
        # np.frombuffer is quicker than unpack
        floatData = np.frombuffer(binData,dtype=np.float32)
        print "Reading dataL ", filename, "block length: ", blockLen, "simData shape", floatData.shape
        binfile.close()
    return floatData.reshape((nsamp,nchan))

def readSetting(filename):
    """
    read settings from yaml file
    """
    # reading telescope setting from simulate file
    with open(filename, 'r') as fp:
        params=yaml.load(fp.read(), Loader=yaml.FullLoader)
    return params

def smooth_bandpass(bandpass):
    """
    get bandpass from fits file
    smooth the bandpass
    """
    output_nos=[]
    output_bandpass=[]
    def smooth(sig,threshold=4.5, level=7, wavelet='db8'):
        sigma = sig.std()
        dwtmatr = pywt.wavedec(data=sig, wavelet=wavelet, level=level)
        denoised = dwtmatr[:]
        denoised[1:] = [pywt.threshold(i, value=threshold*sigma, mode='soft') for i in dwtmatr[1:]]
        smoothed_sig = pywt.waverec(denoised, wavelet, mode='sp1')[:sig.size]
        noises = sig - smoothed_sig
        return smoothed_sig, noises
#    print "bandpass.shape",bandpass.shape 
    idxbad_chan=[]
    inf=float('inf')
#    print "bandpass.shape",bandpass.shape
#    print "length(bandpass)",len(bandpass)
    idxgood_num_=[]
    idxbad_chan_=[]
    sig,nos = smooth(bandpass)
    return sig


def dedisperse(data, dm, chan_freqs, tsamp):
    nt, nf = data.shape
    #delay_time = (4148808.0* dm* (1 / (chan_freqs[-1]) ** 2 - 1 / (chan_freqs) ** 2)/ 1000)
    delay_time = (4148808.0* dm* (1 / (chan_freqs[-1]) ** 2 - 1 / (chan_freqs) ** 2)/ 1000)
    delay_bins = np.round(delay_time / tsamp).astype("int64")
    dedispersed = np.zeros(data.shape, dtype=np.float32)
    for ii in range(nf):
        dedispersed[:, ii] = np.concatenate([data[-delay_bins[ii] :, ii],data[: -delay_bins[ii], ii],])
    return dedispersed
