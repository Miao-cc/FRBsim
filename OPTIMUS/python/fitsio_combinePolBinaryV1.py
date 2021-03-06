import os
import sys
import time
import math
import pywt
import fitsio
import datetime
import numpy as np 
#from pylab import *
from array import array
from decimal import Decimal
import astropy.io.fits as pyfits
import beamWeightMould
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from scipy import interpolate
#import pandas as pd



#smooth the bandpass
#need change#
def smooth_bandpass(bandpass):
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




######################################################################
# 2018/04/04 #  mcc       output: *.fits (pol averaged obs data + binary data)
######################################################################


Flag = 0
Scale = 1
Percent = 0.5

print 'record start time:'
starttime=datetime.datetime.now()
print starttime

#get input files
if (len(sys.argv) == 4):
    infile=sys.argv[1]
    rowdatafile=sys.argv[2]
    outfile=sys.argv[3]
    current_path = os.path.dirname(__file__)
    beamWeightFile = current_path+'/noBeamWeight.dat'

elif (len(sys.argv) == 5):
    infile         = sys.argv[1]
    rowdatafile    = sys.argv[2]
    beamWeightFile = sys.argv[3] 
    outfile        = sys.argv[4]

else:
  print 'too few inputs!'
  print 'example:'
  print 'usage: python fits_combinePolBinary.py fitsfile binaryfile combinefile'
  print 'usage: python fits_combinePolBinary.py fitsfile binaryfile beamWeightfile combinefile'
  sys.exit()





#u19700101=62135683200.0

#==============================================================
#get data in fits file
fits=fitsio.FITS(infile)

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
nchan=header0['OBSNCHAN']
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
print 'NSBLK: ',nsblk,'sample time(s): ',tbin,'channel width(MHz): ',chan_bw

nline=header1['NAXIS2']
pol_type=header1['POL_TYPE']



rmcomm='rm -f '+outfile
os.system(rmcomm)
print rmcomm,outfile
fitsout=fitsio.FITS(outfile,'rw')

dataout = np.zeros(1,dtype=[('TSUBINT','float64'),('OFFS_SUB','float64'),('LST_SUB','float64'),('RA_SUB','float64'),('DEC_SUB','float64'),('GLON_SUB','float64'),('GLAT_SUB','float64'),('FD_ANG','float32'),('POS_ANG','float32'),('PAR_ANG','float32'),('TEL_AZ','float32'),('TEL_ZEN','float32'),('DAT_FREQ','float32',(nchan)),('DAT_WTS','float32',(nchan)),('DAT_OFFS','float32',(nchan)),('DAT_SCL','float32',(nchan)),('DATA','uint8',(nsblk,1,nchan,1))])


fchannel = fits[1].read(rows=[0], columns=['DAT_FREQ'])[0][0][0:nchan]
#print 'DAT_FREQ:', fchannel, fchannel.shape

#==============================================================
#get bandpass and smooth it
#bandpass=get_bandpass(fits)
#smoothBandpass=smooth_bandpass(bandpass)
#smoothBandpass=smoothBandpass/(52.4288/0.0002)
#print "normalized ",'max(smoothBandpass):',np.max(smoothBandpass),'min(smoothBandpass):',np.min(smoothBandpass)
#plot(bandpass,',')
#plot(smoothBandpass,',')
#show()



#==============================================================
# read binanry data and use smoothed bandpass changing the profile in each channel
#binary recorded from float32 to uint8
#simdata=np.zeros((64,nsblk,1,nchan,1))

#print "binarydata.shape:",rowdata.shape
#print "simdata.shape",simdata.shape

#rowdata=np.fromfile(rowdatafile,dtype=np.float32,count=-1)
#for i in range(64):
#    for j in range(nsblk):
#        simdata[i,j,0,:,0]=rowdata[(i*4096+j)*nchan:(i*4096+j+1)*nchan]


print "start reshape file",datetime.datetime.now()
simdata = np.fromfile(rowdatafile,dtype=np.float32,count=-1).reshape((nline,nsblk,1,nchan,1),order='C')
Weight = np.fromfile(beamWeightFile,dtype=np.float64,count=-1).reshape((nline,nsblk,1,nchan,1))
print "end reshape file",datetime.datetime.now()


print 'simdata.dtype',simdata.dtype,'simadata.max',np.max(simdata),'simdata.min',np.min(simdata)

specindex = -1.6
specshap = (fchannel/270.)**(specindex)
specshap[fchannel<270] = 0.
#for i in range(nchan):
    #simdata[:,:,0,i,0]=simdata[:,:,0,i,0]*smoothBandpass[i]*(fchannel[i]/270.)
#simdata[:,:,0,:,0]*=(smoothBandpass*specshap)

print 'simdata.dtype',simdata.dtype,'simadata.max',np.max(simdata),'simdata.min',np.min(simdata)
#simdata[simdata > 127.] = 127.
#simdata[simdata < 0.] = 0.

#simdata = simdata.astype(uint8)
#print 'simdata.dtype',simdata.dtype,'simadata.max',np.max(simdata),'simdata.min',np.min(simdata)

#simdata = np.fromfile(rowdatafile,dtype=np.float32,count=-1).reshape((64,nsblk,1,nchan,1))




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
print 'data.dtype',data[0][0].dtype,'data.max',np.max(data[0][0]),'data.min',np.min(data[0][0])
#==============================================================
#for subint 1 : add the binary data and real obs data
#for subindex in range(nsblk):
    #temp = data[0][0][subindex,0,:,0]/2
    #temp[temp> 127] = 127
    #temp[temp< 0 ] = 0
    #dataout['DATA'][0][subindex,0,:,0] = temp +simdata[rowindex,subindex,0,:,0]

tmpdata = data[0][0][:,0,:,0]
smoothBandpass = smooth_bandpass(tmpdata.sum(axis=0))/nsblk
fidx = np.arange(len(fchannel))
fi = fidx[fchannel >= 270.].min()
smoothBandpass /= smoothBandpass[fi]
tmpdata += (simdata[rowindex,:,0,:,0] * (smoothBandpass*specshap)* Weight[rowindex,:,0,:,0]).astype('uint8') 
#tmpdata /= 2# set max(tmpdata> 255) =255 rather than tmpdata/2
tmpdata[tmpdata > 255] = 255
tmpdata[tmpdata < 0 ] = 0
dataout['DATA'][0][:,0,:,0] = tmpdata

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
    #for subindex in range(nsblk):
        #temp = data[0][0][subindex,0,:,0]/2 
        #temp[temp> 127] = 127
        #temp[temp< 0 ] = 0
        #dataout['DATA'][0][subindex,0,:,0] = temp+simdata[rowindex,subindex,0,:,0]

    tmpdata = data[0][0][:,0,:,0]
    #plot(smoothBandpass)
    #show()
    smoothBandpass = smooth_bandpass(tmpdata.sum(axis=0))/nsblk
    fidx = np.arange(len(fchannel))
    fi = fidx[fchannel >= 270.].min()
    smoothBandpass /= smoothBandpass[fi]
    tmpdata += (simdata[rowindex,:,0,:,0] * (smoothBandpass*specshap)* Weight[rowindex,:,0,:,0]).astype('uint8') 
    #tmpdata /= 2 # set max(tmpdata> 255) =255 rather than tmpdata/2
    tmpdata[tmpdata > 255] = 255
    tmpdata[tmpdata < 0 ] = 0
    dataout['DATA'][0][:,0,:,0] = tmpdata

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
