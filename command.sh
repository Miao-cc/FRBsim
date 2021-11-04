root=`pwd`

# create the folder
cd ${root}
mkdir simBinarydata fastdata

#  create the band spec
python bandspec.py 500 2000
python bandspec.py 500 2000
python bandspec.py 1300 1700
python bandspec.py 1300 1600
python bandspec.py 1400 1600

# simulate the single pulse
cd ${root}/simBinarydata
${root}/OPTIMUS/simulateSimplePsr_mcc -p ${root}/paramsFAST_19Beam.params -p ${root}/params/test_DM_200_W10_0.001.params -o test_DM_200_W10_0.001.dat

# simulate the FRB
cd ${root}/
python simMultiPsr.py simulate.yaml
