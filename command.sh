root=`pwd`

cd ${root}
mkdir simBinarydata fastdata

# 
python bandspec.py 500 2000
python bandspec.py 500 2000
python bandspec.py 1300 1700
python bandspec.py 1300 1600
python bandspec.py 1400 1600

# 
ln -s params/FAST_19Beam.params simBinarydata/
cd ${root}/simBinarydata
${root}/OPTIMUS/simulateSimplePsr_mcc -p FAST_19Beam.params -p ../params/test_DM_200_W10_0.001.params -o test_DM_200_W10_0.001.dat

# 
cd ${root}/
python simMultiPsr.py simulate.yaml
