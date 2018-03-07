#!/usr/bin/env bash
# in run analysis
export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
module use ${MUGQIC_INSTALL_HOME}/modulefiles
module load mugqic/python/2.7.11

python2.7 ../post_pipeline.py "../model/experiments/u0.01875/1_0_0_outs_10/"
python2.7 ../post_pipeline.py "../model/experiments/u0.01875/1_0_005_outs_10/"
python2.7 ../post_pipeline.py "../model/experiments/u0.01875/1_0_01_outs_10/"
python2.7 ../post_pipeline.py "../model/experiments/u0.01875/1_0_02_outs_10/"
python2.7 ../post_pipeline.py "../model/experiments/u0.01875/1_0_065_outs_10/"

