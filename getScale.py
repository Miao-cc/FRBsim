import sys
import time
import fitsio
import datetime
import numpy as np 
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from simFRB import *


print 'record start time:'
starttime=datetime.datetime.now()
print starttime

Num = 1

SettingFile = sys.argv[1]

SettingFile = usage(SettingFile)
# read input parameters
params = readSetting(SettingFile)

for keys in params:
    print "    ", keys, ":", params[keys]

specFile = params["specdata"]
binaryFile = params["binarydata"]
fastFile = params["fastdata"]
outFile = params["outfile"]
scale =  params["scale"]

tele_param = readSetting("params/FAST_19Beam.params")
simBinaryLen = tele_param["t1"]
simBinaryTsamp = tele_param["tsamp"]
simBinaryNchan = tele_param["nchan"]

simNsamp = int(simBinaryLen / simBinaryTsamp)
#simdata=np.zeros(simNsamp*1*simBinaryNchan*1)

print "Simulation File Info: length: %s sec, sample time: %s sec, number of channel: %s" %(simBinaryLen, simBinaryTsamp, simBinaryNchan)

# read header
fits=fitsio.FITS(fastFile[0], mode="r", memmap=True)
hdu0 = fits[0]
hdu1 = fits[1]
header0 = hdu0.read_header()
header1 = hdu1.read_header()
obsfreq=header0['OBSFREQ']
obsbw=header0['OBSBW']
nchan=header0['OBSNCHAN']
nsblk=header1['NSBLK']
tbin=header1['TBIN']
chan_bw=header1['CHAN_BW']
nline=header1['NAXIS2']
npol = header1['NPOL']
fchannel = fits[1].read(rows=[0], columns=['DAT_FREQ'])[0][0][0:nchan]

if simBinaryNchan!=nchan:
    print "Error in channel number"
    sys.exit(0)

# define the read subint nunber
nsubint = int(simNsamp // nsblk)
print nsubint

# load the noise file 
# define the simulate data
noiseData=hdu1.read(rows=range(nsubint), columns=['DATA'])
noiseData = np.array([noiseData[i][0] for i in range(nsubint)])
simdata=np.zeros(nsubint*nsblk*simBinaryNchan*1)

fits.close()

print 'DAT_FREQ:', fchannel, fchannel.shape
print 'NSBLK: ', nsblk, 'sample time(s): ', tbin, 'channel width(MHz): ', chan_bw, "NAXIS2: ", nline
print 'noiseData.dtype',noiseData[0][0].dtype,'noiseData.max',np.max(noiseData[0][0]),'noiseData.min',np.min(noiseData[0][0]),"noiseData.shape", noiseData.shape, "Time sacle (sec): ", tbin*nsubint*nsblk


noiseData = (np.swapaxes(noiseData.squeeze().reshape((-1, npol, nchan)), 1, 2)).sum(axis=2) / 2.
noiseSTD = np.std(noiseData)
noiseMean = np.mean(noiseData)
print noiseMean, noiseSTD


specshape = np.load(specFile[Num])
# create an array only have 4 subint
# add offset in subint 
print "start reshape file",datetime.datetime.now()
# get random offset for diff file
with open(binaryFile[Num], 'rb') as binfile:
    # def read block
    blockLen = nsblk*nchan*(nsubint)
    print blockLen
    # read block
    binData = binfile.read(4*blockLen)
    # np.frombuffer is quicker than unpack
    floatData = np.frombuffer(binData,dtype=np.float32)
    print floatData.shape
    
    # def array start and end, start+offset, end+offset
    arrayStart = 0*blockLen
    arrayEnd = (1+0)*blockLen
    #simdata[arrayStart:arrayEnd] = floatData*scale
    simdata[arrayStart:arrayEnd] = floatData
    
    print "Reading File: ", binaryFile, " block len ", blockLen
    binfile.close()

simdata = simdata.reshape((nsubint*nsblk,nchan))

print "end reshape file",datetime.datetime.now()
print 'simdata.dtype',simdata.dtype,'simadata.max',np.max(simdata),'simdata.min',np.min(simdata)

simdata = simdata*noiseSTD*3

## dedisperse
outdata =(noiseData + np.multiply(simdata, specshape)).astype('uint8')
#data_dedis = dedisperse(outdata, 565, fchannel, tbin)
data_dedis = dedisperse(outdata, 200, fchannel, tbin)
print data_dedis.max()
print data_dedis.shape



# if True, plt the distribution
plot = True
#plot = False
if (plot):
    plt.subplot(3,2,1)
    plt.imshow(simdata.T, origin='lower', aspect='auto', cmap='binary')
    plt.ylim(0,nchan)
    plt.title("Simulated pulse")
    plt.ylabel("Nchan")

    plt.subplot(3,2,3)
    plt.imshow(np.multiply(simdata, specshape).T, origin='lower', aspect='auto', cmap='binary')
    plt.ylim(0,nchan)
    plt.title("Spectrum shaped simulation data")
    plt.ylabel("Nchan")

    plt.subplot(3,2,5)
    plt.imshow(outdata.T, origin='lower', aspect='auto', cmap='binary')
    plt.ylim(0,nchan)
    plt.title("Simulation data + Backend noise")
    plt.ylabel("Nchan")
    plt.xlabel("time sample")

    plt.subplot(3,2,2)
    plt.imshow(noiseData.T, origin='lower', aspect='auto', cmap='binary')
    plt.ylim(0,nchan)
    plt.title("Backend noise")
    plt.ylabel("Nchan")

    plt.subplot(3,2,4)
    plt.imshow(data_dedis.T, origin='lower', aspect='auto', cmap='binary')
    plt.ylim(0,nchan)
    plt.title("Dedisperse simData")
    plt.ylabel("Nchan")

    plt.subplot(3,2,6)
    pulseProfile = data_dedis.sum(axis=1)
    plt.plot(pulseProfile)
    plt.xlim(0,len(pulseProfile))
    plt.title("Pulse profile")
    plt.xlabel("time sample")

    plt.tight_layout()
    plt.savefig("check_plot.png", dpi=300)
    #plt.show()
