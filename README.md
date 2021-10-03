1. create a weight file in freq
 - Usage: start_freq and end_freq in MHz
> python bandspec.py 1150 1250

![bandspec](bandspec.png)

2. create a binary file
3. 
 ```
 cd OPTIMUS
 make
```

3. simulate the file
  > OPTIMUS/simulateSimplePsr_mcc -p parmas/FAST_19Beam.params -p params/test_DM_300_P0_0.0007236s.params -o test_DM_300_P0_0.0007236s.dat

Example:
```
name: J1950+30
p0: 1.610612736
dm: 565
raj: 4.510914803
decj: 0.136026590
flux: 0.1
width: 0.007893
```
This example will simulate a single pulse with pulse width equals to `1.610612736*0.007893=0.012712566325248 sec (P0*width)`

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
    - fastdata:


5. combine the simulation data and backend file

 ` python simMultiPsr.py simulateFRB1000bursts.yaml       (MiaoCC) `
 
 ` python simMultiPsr_RMS.py simulateFRB1000bursts.yaml   (WangPei)`

6. plot for checking the simulating results
`python getScale.py simulate.yaml`
 
 DM=565,,P=1.6s,simBinaryData length = 4 sec

