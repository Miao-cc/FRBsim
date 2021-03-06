# take care
  By now. we fixed the bugs in the code simulateSimplePsr.c
and we write the new code in simulateSimplePsr_mccV1.c
This code can simulate the simple pulsar in each dm and period.
By the way, it used openmp to accelerate the simulate and all of the thread will be used.

# OPTIMUS
1.web url:
  http://www.atnf.csiro.au/people/ghobbs/optimus/index.html?courseSelect=4

2.complie:
  gcc -O3 -lm -o inspectBinaryFile inspectBinaryFile.c simulate.c T2toolkit.c
  gcc -O3 -lm -o simulateSystemNoise simulateSystemNoise.c simulate.c T2toolkit.c
  gcc -O3 -lm -o simulateCal simulateCal.c simulate.c T2toolkit.c
  gcc -lm -o createSearchFile createSearchFile.c simulate.c T2toolkit.c -L{$CFITSIO}/lib -lcfitsio -O3
  gcc -O3 -lm -o simulateSimplePsr simulateSimplePsr.c simulate.c T2toolkit.c
  gcc -O3 -lm -o simulateRFI simulateRFI.c simulate.c T2toolkit.c
  gcc -O3 -lm -o simulateComplexPsr simulateComplexPsr.c simulate.c T2toolkit.c t1polyco.c tempo2pred.c cheby2d.c

  -L{$CFITSIO}/lib is the dir PATH contanting the cfitsio C include file

example:
  gcc -lm -o createSearchFile createSearchFile.c simulate.c T2toolkit.c -L/home/m/psrsoft/cfitsio/lib -lcfitsio -O3

3.Overview

Basic usage of the software

a.Overview of simulating pulsar search mode data files

  Search-mode data files can be simulated using the simulateSearch package. As 
with almost all pulsar software packages, simulateSearch is still unde 
development. The software provides tools for the user to simulate:
 - pulse sequences from pulsars
 - calibration sources
 - receiver and sky noise
 - a particular survey strategy (i.e., pointed observations or drift-scans)
 - different receivers (single pixel, wide-band, multibeams, etc.)
 - backends (specific time sampling, number of frequency channels, digitisation etc.)
  The different sources are simulated independently and then those signals can 
be added together with different simulated backend systems or with different 
observing strategies. The resulting data files are PSRFITS-format search mode 
data files.


b.Simulating system noise

Simulating system noise is trivial:

 optimus]$ ./simulateSystemNoise -o snoise.dat

  This will use default options for the properties of the noise (in this case 
the default is to produce Gaussian random noise with a mean of 0 and a standard 
deviation of 1). 
  The "-o" command line argument provides the output data file (in this case 
snoise.dat). This is a binary data file and cannot be directly viewed. This 
binary file is required for later processing, but it is useful for the user 
simply to check that the simulated noise looks reasonable. To do this we can 
write a small part of the data in text format using the "text output" command 
line argument (-to):

 optimus]$ ./simulateSystemNoise -o snoise.dat -to snoise.txt 0 3 0 3

  This will write, in plain text format, to the file snoise.txt. As the data 
file will become huge we have requested only a time range between 0 and 3 
seconds and frequency channels between 0 and 3 (the last two values on the 
command line). 
  This text file can be plotted in e.g., gnuplot to give:

    gnuplot> plot "snoise.txt" using 1:2

  This should give an output similar to: figure1.pdf

  We confirm that the data is noise like, has a mean of zero and rms of 1. The 
time axis (in seconds) represents the requested three seconds of data.
  Note that neither the binary data file nor the text file can be used in 
standard pulsar software packages such as presto. To do that we need to convert 
the binary data file into PSRFITS format.


c.Keeping track of the binary files
  The simulateSearch software package produces binary data files. These files 
can become confusing if a large number are made. In order to inspect the 
content of a particular file you can use the inspectBinaryFile code:

 optimus]$ ./inspectBinaryFile -f snoise.dat

  The output will list the time span of the data file, the sampling time, the 
observing band, number of channels and whether this particular file corresponds 
to a particular sky direction (note that pulsar signals clearly relate to the 
position of the pulsar, whereas system noise is not position-dependent).
  A typical output is:


            Format:       FORMAT 1
            Name:         Parkes system
            t0 (sec):     0.000000
            t1 (sec):     30.000000
            tsamp (sec):  0.00025
            f1 (MHz):     1518.000000
            f2 (MHz):     1230.000000
            nchan:        96
            RAJ (rad):    0
            DECJ (rad):   0
            Use angle:    0
            Random seed:  -1499991522


d.Creating a PSRFITS search mode file

  In order to use the output of the simulation in standard pulsar search-mode 
visualisation or processing packages we need to produce an output file in the 
PSRFITS search mode format. The PSRFITS search-mode data files are usually 
recorded as 1, 2, 4 or 8 bit values. Using such low bit numbers requires a way 
to define exactly what the individual bits mean. The simulateSearch software 
does not do this in an ideal manner yet. However, for this simple example let's 
make a simple 1-bit data file:

 optimus]$ ./createSearchFile -o search.sf -f snoise.dat 

  This command will load in the binary data file (snoise.dat) and produce a 
PSRFITS search mode file with filename "search.sf".
  The createSearchFile software requires information on the PSRFITS format. 
This comes from a file: psrheader.fits. For now this file (which comes with the
software distribution) needs to be in the local directory when running the code
  Standard visualisation tools can be used to view this PSRFITS search-mode 
file. A simple plot can be produced using

> pfits_plot -f search.sf -s1 1 -s2 1
This should produce a plot similar to: figure2.png

Here the x-axis represents the time samples and the y-axis represents the 
frequency channels. Each point is either a 0 (white) or a 1 (yellow). As only 
system noise was simulated the pulsar search software (presto/sigproc) should 
not be able to find a pulsar in this data set.


e.Adding a simple pulsar

  The software provides two methods for simulating pulsars - a quick and easy 
way (as described here) and a more complex (and slower) method. To create a 
default pulsar we can type:

 optimus]$ ./simulateSimplePsr -o pulsar.dat

  This produces a binary data file (pulsar.dat)
  We can now produce a PSRFITS search mode file containing the system noise and 
the pulsar signal using:

 optimus]$ ./createSearchFile -o search.sf -f snoise.dat -f pulsar.dat

  We can view this with

  > pfits_plot -f search.sf -s1 1 -s2 1

to get a plot similar to the following in which the pulses from the simulated 
pulsar are clearly detectable:

  The default amplitude for the pulses is 3 (in arbitrary units), a period of 
0.3 seconds and a dispersion measure 50cm-3pc.


4.A more realistic simulation

a.Defining the observing system

  In the previous examples we have used the default number of frequency channels
, time sampling, number of bits for recording the data etc. The pulsar we 
simuated had a default period and dispersion measure. All these parameters can 
be changed using parameter files (to save confusion note that this parameter 
file has nothing to do with a ".par" file used in pulsar timing).
  Let's create an initial parameter file with filename parkes1.params with the 
following content:

    name: Parkes system
    f1: 1518
    f2: 1230
    nchan: 128
    t0: 0
    t1: 30
    raj: 0
    decj: 0
    useAngle: 0

  This parameter file has a particular name (Parkes system), a specified 
frequency coverage (from 1518MHz to 1230MHz) and a 128 frequency channels.
  We can now simulate system noise using these input parameters as follows:

 optimus]$ ./simulateSystemNoise -o snoise.dat -p parkes1.params

  Note that we can only combine binary data files (using the createSearchFile 
code) that have the same number of channels, data span etc.

b.Setting the gain and system temperature

  Our simulation does not, so far, include any knowledge of the telescope gain 
nor the system temperature. These parameters can be set in a parameter file 
using:

	gain: (value in K/Jy)
	tsys: (system temperature in K)
    
  For Parkes we can use values similar to:

	gain: 1
	tsys: 25
    
with these parameters set "simulateSystemNoise" will produce noise based on the 
radiometer equation.

c.Setting the number of bits in the PSRFITS output file

  The parameter file can be used to set the number of bits in the output 
PSRFITS file:

	nbits: (number of bits)

  The number of bits can currently be 1, 2 or 8. For example:

        nbits: 2 

  The bits are set as follows:

 -  For 1 bit data, the bit is set to 0 if the value is negative and 1 if the value is positive.
 - For 2 bit data, the levels are currnetly hardcoded to 1= less than -80, 2 = -80 < value < 0, 3 = 0 < value < 80, 4 = value > 80. Note that this needs updating
 - For 8 bit data, the levels (between 0 and 255) are hardcoded to be determined from value*700+100. Note that this needs updating.

  As noted above, the level setting procedure is not yet advanced and will be updated soon

d.Setting the properties of a simple pulsar

  A parameter file can be used to define the properties of a pulsar. An example
is given below (in file pulsar1.params):

	name: Pulsar1
	p0: 0.1
	dm: 20
	raj: 4.510914803
	decj: 0.13602659
	width: 0.03
	flux: 10
	useAngle: 0
    
  The parameters here are:
   - p0: pulse period in seconds
   - dm: dispersion measure in cm-3pc
   - raj: right ascension in radians
   - decj: declination in radians
   - width: pulse width in seconds
   - flux: flux density in Jy
  The pulsar can then be simulated using:

 optimus]$ ./simulateSimplePsr -p parkes1.params -p pulsar1.params -o pulsar.dat
    
  Note that the flux density should only be set if the system noise is being 
correctly determined using the radiometer equation.

e.Simulating two pulsars in a globular cluster

  We wish to simulate two simple pulsars in a globular cluster (i.e., they both 
can be detected in a single, pointed observation). We must first set up out 
observing system. Let's put the following into parkes1.params:

       name: Parkes system
       f1: 1518
       f2: 1230
       nchan: 128
       t0: 0
       t1: 30
       raj: 0
       decj: 0
       useAngle: 0
    
  and now we can simulate the system noise:

 optimus]$ ./simulateSystemNoise -o snoise.dat -p parkes1.params
    
  Now let's simulate the first pulsar with a dispersion measure of 20cm-3pc and 
a pulse period of 0.3 seconds. Let's make the following parameter file (called 
pulsar1.params):

	name: pulsar1
	p0: 0.3
	dm: 20
	raj: 4.510914803
	decj: 0.13602659
	flux: 0.3
	width: 0.05
	useAngle: 0
    
  We can simulate this pulsar using

 optimus]$ ./simulateSimplePsr -p parkes1.params -p pulsar1.params -o pulsar1.dat
    
  Now for the second pulsar (in pulsar2.params):

        name: pulsar2
        p0: 0.08
        dm: 50
        raj: 4.510914803
        decj: 0.13602659
        flux: 1
        width: 0.01	
        useAngle: 0
    
  We can simulate this pulsar using

 optimus]$ ./simulateSimplePsr -p parkes1.params -p pulsar2.params -o pulsar2.dat
    
  We should now have a binary data file for the system noise (snoise.dat) and 
two binary data files for the two pulsars (pulsar1.dat and pulsar2.dat). 
Finally we have to make the output PSRFITS file. Let's assume that we want to 
write the data as 1-bit values. To do this let's make another parameter file 
(digitiser.params) with:

         name: digitiser1
         nbits: 1
   
  and now we can create the FITS file:

 optimus]$ ./createSearchFile -o search.sf -p parkes1.params -p digitiser.params -f snoise.dat -f pulsar1.dat -f pulsar2.dat   

  That should produce the PSRFITS file (search.sf) that can be viewed using:

> pfits_plot -f search.sf -s1 1 -s2 1
     
 This should produce a result similar to:

  The thick lines correspond to the pulses from pulsar 1. The thin lines from 
pulsar 2.

5.Using sky positions
a.Introduction to using sky positions with the simulateSearch software 
  Some signals to be simulated do not depend upon where the virtual telescope 
is pointing (for instance, the receiver noise). However, other signals depend 
strongly on the sky position. For instance, a pulsar will be most detectable 
when the telescope is pointing straight at the pulsar. The signal from the 
pulsar if the telescope is not pointing directly at the pulsar will depend upon 
the shape of the "beam". If the source and beam position should be taken into 
account during the simulation then the parameter file should include:


         useAngle: 1
   
  Conversely if it should not be used then the parameter file should include:

         useAngle: 0
   
  When a source is being simulated and the position must be accounted for then 
its right ascension and declination must be provided (in radians):

         raj: (source right ascension in radians)
         decj: (source declination in radians)
   
  Of course we also need to define where the telescope beam is (or the 
telescope beams are) pointing. We can do this by adding the beamRA0 and 
beamDEC0 parameters into the file used when running createSearchFile:

         beamRA0: (beam right ascension in radians)
         beamDEC0: (beam declination in radiation)
   
  These correspond to the beam pointing position at the start of the observation
. By default the telescope is assumed to track this position and so the beam 
position (in right ascension and declination) will not change.



  - simulateCal???

  - simulateRFI???


  - simulateSimplePsr???
    function: (usage1) simulate the pulsar with the defaultparameters.
                       This produces a binary data file (pulsar.dat)
              (usage2) simulate the pulsar with the parameters in the parameter file (called 
                       pulsar1.params)
    usage1: ./simulateSimplePsr -o [binary file *.dat]
            ./simulateSimplePsr -o pulsar.dat
    usage2: ./simulateSimplePsr -p [observing sys parameter file *.params] -p [pulsar parameter file *.params] -o [binary file *.dat]
            ./simulateSimplePsr -p parkes1.params -p pulsar2.params -o pulsar2.dat
    out: *.dat

  - simulateComplexpsr???


  - simulateSystemNoise???
    function: (usage1) will use default options for the properties of the noise (in this case 
              the default is to produce Gaussian random noise with a mean of 0 and a standard 
              deviation of 1).
              (usage2) This will write, in plain text format, to the file snoise.txt.
              (usage3) simulate system noise using input parameters in the *.params files
    usage1: ./simulateSystemNoise -o [binary file *.dat]
            ./simulateSystemNoise -o snoise.dat
    usage2: ./simulateSystemNoise -o [binary file *.dat] -to [ASIC file *.txt] [t0] [t1] [chan0] [chan1]
            [t0]: [start fo time]; [t1]: [end fo time]; ===> time range between 0 and 3 seconds 
            [chan0]: [start num fo chnnels]; [chan1]: [end um fo chnnels]; ===> frequency channels between 0 and 3
            ./simulateSystemNoise -o snoise.dat -to snoise.txt 0 3 0 3
    usage3: ./simulateSystemNoise -o [binary file *.dat] -p [parameter file *.params]
            ./simulateSystemNoise -o snoise.dat -p parkes1.params
    out: *.txt ; *.dat

  - createSearchFile???
    function: will load in the binary data file (snoise.dat) and produce a 
              PSRFITS search mode file with filename "search.sf".
    usage: ./createSearchFile -o [outfile *.sf] -f [binary file *.dat]       
           ./createSearchFile -o search.sf -f snoise.dat 
           ./createSearchFile -o [outfile *.sf] -p [parameter file *.params] -f [binary file *.dat]       
           ./createSearchFile -o search.sf -p parkes1.params -p digitiser.params -f snoise.dat -f pulsar1.dat -f pulsar2.dat 
    out: *.sf

  - inspectBinaryFile??? 
    function: view binary files 
    usage: ./inspectBinaryFile -f [binary file *.dat]
           ./inspectBinaryFile -f snoise.dat
    out: will list the time span of the data file, the sampling time, the observing 
         band, number of channels and so on.















