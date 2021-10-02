root=/home/wangpei/FRBsim
cd ${root}

# frb 121102 sim

#cd ${root}/simBinarydata
#${root}/OPTIMUS/simulateSimplePsr_mcc -p FAST_19Beam.params -p test_DM_300_P0_0.005s.params -o test_DM_300_P0_0.005s.dat

#${root}/OPTIMUS/simulateSimplePsr_mcc -p FAST_19Beam.params -p test_DM_200_P0_0.0007236s.params -o test_DM_200_P0_0.0007236s.dat

#cd ${root}/
#python simMultiPsr.py simulate-Test.yaml
#python simMultiPsr.py simulate-Test1.yaml


python bandspec.py 1300 1600
python bandspec.py 1300 1700
python bandspec.py 1300 1500
python bandspec.py 1400 1600
python bandspec.py 1400 1700
python bandspec.py 1200 1600
python bandspec.py 1200 1300
