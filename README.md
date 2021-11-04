1. create a weight file in freq
 - Usage: start_freq and end_freq in MHz
> python bandspec.py 1150 1250

![bandspec](bandspec.png)

2. create a binary file

 ```
 cd OPTIMUS
 make
```

3. simulate the file
  > OPTIMUS/simulateSimplePsr_mcc -p params/FAST_19Beam.params -p params/test_DM_200_W10_0.001.params -o test_DM_200_W10_0.001.dat

Example:
test_DM_200_W10_0.001.params
```
name: J1950+30
p0: 0.8
dm: 565
raj: 4.510914803
decj: 0.136026590
flux: 0.1
width: 0.007893
```
FAST_19Beam.params
```
name: FAST system
f1: 1000
f2: 1500
nchan: 4096
t0: 0.00000
t1: 0.6
gain: 15.00
tsys: 30.00
raj: 0.000000000
decj: 0.000000000
useAngle: 0
tsamp: 0.000049152
digitiser: 8
```

This example will simulate a single pulse with pulse width equals to `0.8*0.007893=0.0063144 sec (P0*width)`

 - Tips: 
  1. pulse start time delay is `pow(2,10)*tsamp = 0.0503316` sec 
  2. `t1` in `FAST_19Beam.params` should be more than time delay cause by DM, like t1=0.6 > (0.461111+0.0503316)
  3. `P0` in file `test_DM_565_W10_0.0078.params` should be more than the `t1` in file `FAST_19Beam.params`, like p1=0.8 > t1=0.6

4. define the pulses

```
fastdata:
  - fastdata/M92_tracking-M01_0214.fits
binarydata:
  - simBinarydata/test_DM_200_P0_0.0007236s.dat
  - simBinarydata/test_DM_200_P0_0.0007236s.dat
  - simBinarydata/test_DM_200_P0_0.0007236s.dat
  - simBinarydata/test_DM_200_P0_0.0007236s.dat
  - simBinarydata/test_DM_200_P0_0.0007236s.dat
specdata:
  - bandspec_1200-1600MHz.npy
  - bandspec_500-2000MHz.npy
  - bandspec_1300-1700MHz.npy
  - bandspec_1300-1600MHz.npy
  - bandspec_1400-1600MHz.npy
toafile:
  - FRB121102TOA.tim
outfile:
  - simulateTest1.fits
scale:
  - 15.0
  - 8.0
  - 7.0
  - 6.0
  - 5.0
```
  * fastdata: the telescope noise
  * binarydata: single pulse simulation file with different pulse width
  * specdata: scale used to scale the pulse in frequency
  * toafile: time of arrival
  * outfile: the output result
  * scale: used to addjust pulse S/N


5. combine the simulation data and backend file

 ` python simMultiPsr.py simulate.yaml       (MiaoCC) `
 
 ` python simMultiPsr_RMS.py simulate.yaml   (WangPei)`

6. plot for checking the simulating results
`python getScale.py simulate.yaml`

### envirement
python2.7

```
astropy==2.0.16
fitsio==1.0.1
matplotlib==2.2.5
numpy==1.16.2
PyWavelets==1.0.3
PyYAML==5.4.1
scipy==1.2.3

```
