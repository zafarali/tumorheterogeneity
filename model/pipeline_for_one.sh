#!/usr/bin/env bash
# in run analysis
export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
module use ${MUGQIC_INSTALL_HOME}/modulefiles
module load mugqic/python/2.7.11

python2.7 ../post_pipeline.py ${1}
python2 ./TEMP_create_mixing_analysis.py ${1}
