#!/usr/bin/env bash

export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
module use ${MUGQIC_INSTALL_HOME}/modulefiles
module load mugqic/python/2.7.11

cd /lb/project/gravel/zafarali_projects/tumorheterogeneity/figures

PATH_TO_ALL_SIMS="../model/experiments/u0.01875/"
SELECTED_TUMOR="10" # seed of the selected tumor
python2 do_power_analysis.py $PATH_TO_ALL_SIMS $SELECTED_TUMOR
python2 create_AFS_figures.py
python2 create_number_of_divisions_table.py > S1table.txt
python2 create_Splots_figures.py $SELECTED_TUMOR
python2 create_cluster_advantage_plots.py
python2 create_fanning_plots.py
python2 create_mixing_analysis.py "$PATH_TO_ALL_SIMS/1_0_" $SELECTED_TUMOR
python2 create_trees.py $PATH_TO_ALL_SIMS $SELECTED_TUMOR

