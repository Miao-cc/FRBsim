import os
import sys
import time
import fitsio
import astropy.io.fits as pyfits
import datetime
import numpy as np 
import matplotlib.pyplot as plt
from simFRB import *

def readBackend(fo, subintStart, subintEnd):
    backend = []
    for rowindex in range(subintStart, subintEnd):
        data=fits[1].read(rows=[rowindex], columns=['DATA'])
        backend.append(data[0][0][:,:,:,:])
    return np.array(backend)


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
fits=fitsio.FITS(fastFile[0], "r")
# open the new output data
fitsout=fitsio.FITS(outFile[0],"rw")

hdu0 = fits[0]
header0 = hdu0.read_header()
hdu1 = fits[1]
header1 = hdu1.read_header()

hdrver=header0['HDRVER']
date=header0['DATE']
ant_x=header0['ANT_X']
ant_y=header0['ANT_Y']
ant_z=header0['ANT_Z']
obsfreq=header0['OBSFREQ']
obsbw=header0['OBSBW']
nchan=int(header0['OBSNCHAN'])
ra=header0['RA']
dec=header0['DEC']
bmaj=header0['BMAJ']
bmin=header0['BMIN']
date_obs=header0['DATE-OBS']
stt_imjd=header0['STT_IMJD']
stt_smjd=header0['STT_SMJD']
stt_offs=header0['STT_OFFS']
stt_lst=header0['STT_LST']
nsuboffs=header1['NSUBOFFS']
nchnoffs=header1['NCHNOFFS']
nsblk=header1['NSBLK']
npol=header1['NPOL']
tbin=header1['TBIN']
chan_bw=header1['CHAN_BW']
nline=header1['NAXIS2']
pol_type=header1['POL_TYPE']
fchannel = fits[1].read(rows=[0], columns=['DAT_FREQ'])[0][0][0:nchan]
print "nline: ", nline, "NSBLK: ", nsblk, "sample time(s): ", tbin, "channel width(MHz): ", chan_bw
print 'DAT_FREQ:', fchannel, fchannel.shape

# define the out put file
dataout = np.zeros(1,dtype=[('TSUBINT','float64'),('OFFS_SUB','float64'),('LST_SUB','float64'),('RA_SUB','float64'),('DEC_SUB','float64'),('GLON_SUB','float64'),('GLAT_SUB','float64'),('FD_ANG','float32'),('POS_ANG','float32'),('PAR_ANG','float32'),('TEL_AZ','float32'),('TEL_ZEN','float32'),('DAT_FREQ','float32',(nchan)),('DAT_WTS','float32',(nchan)),('DAT_OFFS','float32',(nchan)),('DAT_SCL','float32',(nchan)),('DATA','uint8',(nsblk,1,nchan,1))])

# define the simdata to add the pulse together
simdata=np.zeros((nline*nsblk,nchan))
# read the backend data to get the STD

# add offset in subint 
print "start reshape file",datetime.datetime.now()

for num, filename in enumerate(binaryFile):
    #### read the spec index from file
    specshape = np.load(specFile[num])
    SNR = scale[num]
    binData = readBinary(filename, simBinaryNchan, simNsamp)

    # start subint and end subint
    arrayStart = int(pulseOffset[num]/tbin)
    arrayEnd = arrayStart + simNsamp
    #"""
    subintStart = int(arrayStart/nsblk)
    subintEnd = int(arrayEnd/nsblk)
    if subintEnd>nline: subintEnd=nline
    backend = readBackend(fits, subintStart, subintEnd)
    backend_pol0 = backend[:,:,0,:,0]
    scaleSTD = np.std(backend_pol0.mean(axis=2))
    scaleMean = np.mean(backend_pol0)
    #"""
    # simulate pulse * scale * frequence shape
    #simdata[arrayStart:arrayEnd, :] = binData*SNR*specshape + simdata[arrayStart:arrayEnd, :]
    scaledSimdata = binData*SNR*scaleSTD*specshape
    simdata[arrayStart:arrayEnd, :] = scaledSimdata + simdata[arrayStart:arrayEnd, :]

    print "scaled pulse: max, min: ", np.max(scaledSimdata), np.min(scaledSimdata) , "SNR*scaleSTD", SNR*scaleSTD
    print "Reading spectra shape: ", specFile[num], "SNR Value: ", SNR
    print "Reading simulation data: ", num, filename, "Pulse TOA: ", pulseOffset[num]
    print 'backend_pol0.shape', backend_pol0.shape, 'backend_pol0.min, backend_pol0.max',np.min(backend_pol0), np.max(backend_pol0), 'backend_pol0.std,backend_pol0.mean in time serise', scaleSTD, scaleMean
    print "+"*50


simdata = simdata.reshape((nline,nsblk,1,nchan,1))
print 'simdata.dtype',simdata.dtype,'simadata.max',np.max(simdata),'simdata.min',np.min(simdata)
print "end reshape file",datetime.datetime.now()

### get header
#=======================================================================
rowindex=0
dataout['TSUBINT'][0]=fits[1].read(rows=[rowindex], columns=['TSUBINT'])[0][0]
dataout['OFFS_SUB'][0]=fits[1].read(rows=[rowindex], columns=['OFFS_SUB'])[0][0]
dataout['LST_SUB'][0]=fits[1].read(rows=[rowindex], columns=['LST_SUB'])[0][0]
dataout['RA_SUB'][0]=fits[1].read(rows=[rowindex], columns=['RA_SUB'])[0][0]
dataout['DEC_SUB'][0]=fits[1].read(rows=[rowindex], columns=['DEC_SUB'])[0][0]
dataout['GLON_SUB'][0]=fits[1].read(rows=[rowindex], columns=['GLON_SUB'])[0][0]
dataout['GLAT_SUB'][0]=fits[1].read(rows=[rowindex], columns=['GLAT_SUB'])[0][0]
dataout['FD_ANG'][0]=fits[1].read(rows=[rowindex], columns=['FD_ANG'])[0][0]
dataout['POS_ANG'][0]=fits[1].read(rows=[rowindex], columns=['POS_ANG'])[0][0]
dataout['PAR_ANG'][0]=fits[1].read(rows=[rowindex], columns=['PAR_ANG'])[0][0]
dataout['TEL_AZ'][0]=fits[1].read(rows=[rowindex], columns=['TEL_AZ'])[0][0]
dataout['TEL_ZEN'][0]=fits[1].read(rows=[rowindex], columns=['TEL_ZEN'])[0][0]
dataout['DAT_FREQ'][0]=fits[1].read(rows=[rowindex], columns=['DAT_FREQ'])[0][0][0:nchan]
dataout['DAT_WTS'][0]=fits[1].read(rows=[rowindex], columns=['DAT_WTS'])[0][0][0:nchan]
dataout['DAT_OFFS'][0][0:nchan]=fits[1].read(rows=[rowindex], columns=['DAT_OFFS'])[0][0][0:nchan]
dataout['DAT_SCL'][0][0:nchan]=fits[1].read(rows=[rowindex], columns=['DAT_SCL'])[0][0][0:nchan]

data=fits[1].read(rows=[rowindex], columns=['DATA'])
print 'data.dtype',data[0][0].dtype,'data.max',np.max(data[0][0]),'data.min',np.min(data[0][0]), 'data.shape', data[0][0].shape

# nsblk, pol, channel, 1
tmpdata = data[0][0][:,0,:,0]
#########
#smoothBandpass = smooth_bandpass(tmpdata.sum(axis=0))/nsblk
#fidx = np.arange(len(fchannel))
#fi = fidx[fchannel >= 270.].min()
#smoothBandpass /= smoothBandpass[fi]
#smoothBandpass[0:100] = np.mean(smoothBandpass[100:200])
#tmpdata += (simdata[rowindex,:,0,:,0] * (smoothBandpass).astype('uint8')
#########
# remove smooth bandpass
tmpdata = tmpdata.astype('float')
tmpdata += (simdata[rowindex,:,0,:,0])
#tmpdata /= 2# set max(tmpdata> 255) =255 rather than tmpdata/2
tmpdata[tmpdata > 255] = 255
tmpdata[tmpdata < 0 ] = 0
dataout['DATA'][0][:,0,:,0] = tmpdata.astype('uint8')

fitsout.write(dataout)
#=======================================================================

fitsout[0].write_key('HDRVER',hdrver,comment="")
fitsout[0].write_key('FITSTYPE','PSRFITS',comment="FITS definition ")
fitsout[0].write_key('DATE',date,comment="")
fitsout[0].write_key('OBSERVER','FAST_TEAM',comment="Observer name")
fitsout[0].write_key('PROJID','Drift',comment="Project name")
fitsout[0].write_key('TELESCOP','FAST',comment="Telescope name")
fitsout[0].write_key('ANT_X',ant_x,comment="")
fitsout[0].write_key('ANT_Y',ant_y,comment="")
fitsout[0].write_key('ANT_Z',ant_z,comment="")
fitsout[0].write_key('FRONTEND','WIDEBAND',comment="Frontend ID")
fitsout[0].write_key('NRCVR',1,comment="")
fitsout[0].write_key('FD_POLN','LIN',comment="LIN or CIRC")
fitsout[0].write_key('FD_HAND',1,comment="")
fitsout[0].write_key('FD_SANG',0.,comment="")
fitsout[0].write_key('FD_XYPH',0.,comment="")
fitsout[0].write_key('BACKEND','ROACH',comment="Backend ID")
fitsout[0].write_key('BECONFIG','N/A',comment="")
fitsout[0].write_key('BE_PHASE',1,comment="")
fitsout[0].write_key('BE_DCC',0,comment="")
fitsout[0].write_key('BE_DELAY',0.,comment="")
fitsout[0].write_key('TCYCLE',0.,comment="")
fitsout[0].write_key('OBS_MODE','SEARCH',comment="(PSR, CAL, SEARCH)")
fitsout[0].write_key('DATE-OBS',date_obs,comment="Date of observation")
fitsout[0].write_key('OBSFREQ',obsfreq,comment="[MHz] Bandfrequency")
fitsout[0].write_key('OBSBW',obsbw,comment="[MHz] Bandwidth")
fitsout[0].write_key('OBSNCHAN',nchan,comment="Number of channels")
fitsout[0].write_key('CHAN_DM',0.,comment="")
fitsout[0].write_key('SRC_NAME','Drift',comment="Source or scan ID")
fitsout[0].write_key('COORD_MD','J2000',comment="")
fitsout[0].write_key('EQUINOX',2000.,comment="")

fitsout[0].write_key('RA',ra,comment="")
fitsout[0].write_key('DEC',dec,comment="")
fitsout[0].write_key('BMAJ',bmaj,comment="[deg] Beam major axis length")
fitsout[0].write_key('BMIN',bmin,comment="[deg] Beam minor axis length")
fitsout[0].write_key('BPA',0.,comment="[deg] Beam position angle")
fitsout[0].write_key('STT_CRD1','00:00:00.00',comment="")
fitsout[0].write_key('STT_CRD2','00:00:00.00',comment="")
fitsout[0].write_key('TRK_MODE','TRACK',comment="")
fitsout[0].write_key('STP_CRD1','00:00:00.00',comment="")
fitsout[0].write_key('STP_CRD2','00:00:00.00',comment="")
fitsout[0].write_key('SCANLEN',0.,comment="")
fitsout[0].write_key('FD_MODE','FA',comment="")
fitsout[0].write_key('FA_REQ',0.,comment="")
fitsout[0].write_key('CAL_MODE','OFF',comment="")
fitsout[0].write_key('CAL_FREQ',0.,comment="")
fitsout[0].write_key('CAL_DCYC',0.,comment="")
fitsout[0].write_key('CAL_PHS',0.,comment="")
fitsout[0].write_key('STT_IMJD',stt_imjd,comment="Start MJD (UTC days) (J - long integer)")
fitsout[0].write_key('STT_SMJD',stt_smjd,comment="[s] Start time (sec past UTC 00h) (J)")
fitsout[0].write_key('STT_OFFS',stt_offs,comment="[s] Start time offset (D)")
fitsout[0].write_key('STT_LST',stt_lst,comment="[s] Start LST (D)")

fitsout[1].write_key('INT_TYPE','TIME',comment="Time axis (TIME, BINPHSPERI, BINLNGASC, etc)")
fitsout[1].write_key('INT_UNIT','SEC',comment="Unit of time axis (SEC, PHS (0-1),DEG)")
fitsout[1].write_key('SCALE','FluxDen',comment="")
fitsout[1].write_key('NPOL',1,comment="Nr of polarisations")
fitsout[1].write_key('POL_TYPE','AABB',comment="Polarisation identifier")
fitsout[1].write_key('TBIN',tbin,comment="[s] Time per bin or sample")
fitsout[1].write_key('NBIN',1,comment="")
fitsout[1].write_key('NBIN_PRD',0,comment="Nr of bins/pulse period (for gated data)")
fitsout[1].write_key('PHS_OFFS',0.0,comment="Phase offset of bin 0 for gated data")
fitsout[1].write_key('NBITS',8,comment="Nr of bits/datum ")
fitsout[1].write_key('NSUBOFFS',nsuboffs,comment="Subint offset ")
fitsout[1].write_key('NCHNOFFS',nchnoffs,comment="Channel/sub-band offset for split files")
fitsout[1].write_key('NCHAN',nchan,comment="Number of channels")
fitsout[1].write_key('CHAN_BW',chan_bw,comment="[MHz] Channel/sub-band width")
fitsout[1].write_key('NSBLK',nsblk,comment="Samples/row ")
fitsout[1].write_key('EXTNAME','SUBINT  ',comment="name of this binary table extension")
fitsout[1].write_key('EXTVER',1,comment="")


#==============================================================
#for subint 2-64 : add the binary data and real obs data
for rowindex in range(1,nline):
    dataout['TSUBINT'][0]=fits[1].read(rows=[rowindex], columns=['TSUBINT'])[0][0]
    dataout['OFFS_SUB'][0]=fits[1].read(rows=[rowindex], columns=['OFFS_SUB'])[0][0]
    dataout['LST_SUB'][0]=fits[1].read(rows=[rowindex], columns=['LST_SUB'])[0][0]
    dataout['RA_SUB'][0]=fits[1].read(rows=[rowindex], columns=['RA_SUB'])[0][0]
    dataout['DEC_SUB'][0]=fits[1].read(rows=[rowindex], columns=['DEC_SUB'])[0][0]
    dataout['GLON_SUB'][0]=fits[1].read(rows=[rowindex], columns=['GLON_SUB'])[0][0]
    dataout['GLAT_SUB'][0]=fits[1].read(rows=[rowindex], columns=['GLAT_SUB'])[0][0]
    dataout['FD_ANG'][0]=fits[1].read(rows=[rowindex], columns=['FD_ANG'])[0][0]
    dataout['POS_ANG'][0]=fits[1].read(rows=[rowindex], columns=['POS_ANG'])[0][0]
    dataout['PAR_ANG'][0]=fits[1].read(rows=[rowindex], columns=['PAR_ANG'])[0][0]
    dataout['TEL_AZ'][0]=fits[1].read(rows=[rowindex], columns=['TEL_AZ'])[0][0]
    dataout['TEL_ZEN'][0]=fits[1].read(rows=[rowindex], columns=['TEL_ZEN'])[0][0]
    dataout['DAT_FREQ'][0]=fits[1].read(rows=[rowindex], columns=['DAT_FREQ'])[0][0][0:nchan]
    dataout['DAT_WTS'][0]=fits[1].read(rows=[rowindex], columns=['DAT_WTS'])[0][0][0:nchan]
    dataout['DAT_OFFS'][0][0:nchan]=fits[1].read(rows=[rowindex], columns=['DAT_OFFS'])[0][0][0:nchan]
    dataout['DAT_SCL'][0][0:nchan]=fits[1].read(rows=[rowindex], columns=['DAT_SCL'])[0][0][0:nchan]

    data=fits[1].read(rows=[rowindex], columns=['DATA'])

    tmpdata = data[0][0][:,0,:,0]
    #########
    #smoothBandpass = smooth_bandpass(tmpdata.sum(axis=0))/nsblk
    #fidx = np.arange(len(fchannel))
    #fi = fidx[fchannel >= 270.].min()
    #smoothBandpass /= smoothBandpass[fi]
    #########
    # remove smooth bandpass
    #tmpdata += (simdata[rowindex,:,0,:,0] * (smoothBandpass).astype('uint8')
    tmpdata = tmpdata.astype('float')
    tmpdata += (simdata[rowindex,:,0,:,0])
    tmpdata[tmpdata > 255] = 255
    tmpdata[tmpdata < 0 ] = 0
    dataout['DATA'][0][:,0,:,0] = tmpdata.astype('uint8')

    fitsout[-1].append(dataout)

#fitsout.write(dataout)
#==============================================================


fitsout.close()
print '--------------------------------------------'
print '             Finished!                      '

endtime=datetime.datetime.now()
print 'START:',starttime
print 'END:',endtime
duration=endtime-starttime
print 'DURATION:',duration.seconds,' sec'
