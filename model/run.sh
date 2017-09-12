export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
module use ${MUGQIC_INSTALL_HOME}/modulefiles
module load mugqic/python/2.7.11

./${1}.exe ${1}_outs_${2} 1 ${2} ${@:3}
cp specs.json ${1}_outs_${2}/specs.json
python2.7 ../pipeline ${1}_outs_${2}