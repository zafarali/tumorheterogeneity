#!/usr/bin/env bash

export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
module use ${MUGQIC_INSTALL_HOME}/modulefiles
module load mugqic/python/2.7.11

cd /lb/project/gravel/zafarali_projects/tumorheterogeneity/figures

PATH_TO_ALL_SIMS="../model/experiments/u0.01875"

python2 finescaleAFS.py turnover close $PATH_TO_ALL_SIMS
python2 finescaleAFS.py turnover notclose $PATH_TO_ALL_SIMS
python2 finescaleAFS.py noturnover close $PATH_TO_ALL_SIMS
python2 finescaleAFS.py notturnover notclose $PATH_TO_ALL_SIMS
python2 do_power_analysis.py
python2 create_AFS_figures.py
python2 create_number_of_divisions_table.py
python2 create_SvsN_figures.py
python2 create_Splots_figures.py
python2 create_cluster_advantage_plots.py


