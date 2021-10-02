1. 
# create a weight file in freq
# input start_freq and end_freq in MHz
python bandspec.py 1150 1250

2.
# create a binary file
${root}/OPTIMUS/simulateSimplePsr_mcc -p FAST_19Beam.params -p test_DM_300_P0_0.0007236s.params -o test_DM_300_P0_0.0007236s.dat

3.
# combine the simulation data and backend file
python simMultiPsr_RMS.py simulateFRB1000bursts.yaml

4.
# plot for checking
DM=565
P=1.6s
simBinaryData length = 4 sec

python getScale.py simulate-Test.yaml
