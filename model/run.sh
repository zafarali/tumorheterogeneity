## argument 1 == executable name
## argument 2 == mutation rate
## argument 3 == random seed

export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
module use ${MUGQIC_INSTALL_HOME}/modulefiles
module load mugqic/python/2.7.11

./experiments/${2}/${1}.exe "${1}_outs_${3}_TEMP" 1 ${3}
cp specs.json ${1}_outs_${3}_TEMP/specs.json
python2 ../pipeline ${1}_outs_${3}_TEMP
mkdir ./experiments/u0.01/${1}_outs_${3}
mv ${1}_outs_${3}_TEMP/* ./experiments/u0.01/${1}_outs_${3}
cp -r ${1}_outs_${3}_TEMP/*pipe_out* ./experiments/u0.01/${1}_outs_${3}
rm -rf ${1}_outs_${3}_TEMP