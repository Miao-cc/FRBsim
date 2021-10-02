import os
import sys
import time
import copy
#import fitsio
import astropy.io.fits as pyfits
import datetime
import numpy as np 
import matplotlib.pyplot as plt
from simFRB import *



######################################################################
# 2020/11/22 #  mcc       output: *.fits (pol averaged obs data + binary data)
######################################################################
starttime = datetime.datetime.now()
print 'record start time:', starttime

# reading input settings
SettingFile = sys.argv[1]
# check the input
SettingFile = usage(SettingFile)
# read input parameters
params = readSetting(SettingFile)

# print the parameters
for keys in params:
    print "    ", keys, ":", params[keys]

# define 
specFile = params["specdata"]
binaryFile = params["binarydata"]
fastFile = params["fastdata"]
outFile = params["outfile"]
toaFile = params["toafile"]
scale =  params["scale"]

# readting simulate settings from FAST_19Beam.params
tele_param = readSetting("simBinarydata/FAST_19Beam.params")
simBinaryLen = tele_param["t1"]
simBinaryTsamp = tele_param["tsamp"]
simBinaryNchan = tele_param["nchan"]
simNsamp = int(simBinaryLen / simBinaryTsamp)

print "Simulation File Info: length: %s sec, sample time: %s sec, number of channel: %s" %(simBinaryLen, simBinaryTsamp, simBinaryNchan)

########################################
# get pulse TOA, should between 0~nline*nsblk*tbin
print toaFile
pulseOffset = np.loadtxt(toaFile[0])
print "Reading TOAs from file: %s\n TOA number is : %s" %(toaFile, len(pulseOffset))

# check input
PulseNum = len(pulseOffset)
inputCheck(specFile, PulseNum)
inputCheck(binaryFile, PulseNum)
inputCheck(scale, PulseNum)


########################################
# combine settings

# remove the same name simulation file
rmcomm='rm -f '+outFile[0]
os.system(rmcomm)
print "Remove the old data: ", rmcomm

# open backend noise file
#fits=fitsio.FITS(fastFile[0], "r")
noiseFits= pyfits.open(fastFile[0], "readonly", memmap=True)
#outFits  = pyfits.open(outFile[0], "append", memmap=True)

hdu0 = noiseFits[0]
hdu1 = noiseFits[1]

date=header0['DATE']
nchan=header0['OBSNCHAN']
nsblk=header1['NSBLK']
npol=header1['NPOL']
tbin=header1['TBIN']
chan_bw=header1['CHAN_BW']
nline=header1['NAXIS2']
pol_type=header1['POL_TYPE']
hdu1_data = hdu1.data
data = hdu1_data['DATA']
fchannel = hdu1_data[0]['DAT_FREQ']
print "nline: ", nline, "NSBLK: ", nsblk, "sample time(s): ", tbin, "channel width(MHz): ", chan_bw
print 'DAT_FREQ:', fchannel, fchannel.shape

# define the simdata to add the pulse together
simdata=np.zeros((nline*nsblk,nchan))

# add offset in subint 
print "start reshape file",datetime.datetime.now()

for num, filename in enumerate(binaryFile):
    #### read the spec index from file
    specshape = np.load(specFile[num])
    scaleValue = scale[num]
    binData = readBinary(filename, simBinaryNchan, simNsamp)
    print "+"*50
    print "Reading spectra shape: ", specFile[num], "Scale Value: ", scaleValue
    print "Reading simulation data: ", num, filename, "Pulse TOA: ", pulseOffset[num]
    arrayStart = int(pulseOffset[num]/tbin)
    arrayEnd = arrayStart + simNsamp
    # simulate pulse * scale * frequence shape
    simdata[arrayStart:arrayEnd, :] = binData*scaleValue*specshape + simdata[arrayStart:arrayEnd, :]
simdata = simdata.reshape((nline,nsblk,1,nchan,1))
print 'simdata.dtype',simdata.dtype,'simadata.max',np.max(simdata),'simdata.min',np.min(simdata)
print "end reshape file",datetime.datetime.now()

#=======================================================================
#==============================================================
#add the binary data and real obs data
# sum data
#data = np.expand_dims(np.sum(data[:, :, :, :, 0], axis=2,keepdims=True)/2., axis=5)
#data = np.expand_dims(np.sum(data[:, :, :, :, 0], axis=2,keepdims=True)/2., axis=5)

for rowindex in range(0,nline):
    tmpdata = data[rowindex, :, :, :, :] + simdata[rowindex, :, :, :, :]
    tmpdata += (simdata[rowindex, :, :, :, :]).astype('uint8')
    tmpdata[tmpdata > 255] = 255
    tmpdata[tmpdata < 0 ] = 0
    data[rowindex, :, :, :, :] = tmpdata

noiseFits.writeto('newtable.fits')

#fitsout.write(dataout)
#==============================================================


noiseFits.close()
print '--------------------------------------------'
print '             Finished!                      '

endtime=datetime.datetime.now()
print 'START:',starttime
print 'END:',endtime
duration=endtime-starttime
print 'DURATION:',duration.seconds,' sec'
