#!/usr/bin/env bash

export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
module use ${MUGQIC_INSTALL_HOME}/modulefiles
module load mugqic/python/2.7.11

cd /lb/project/gravel/zafarali_projects/tumorheterogeneity/figures

PATH_TO_ALL_SIMS="../model/experiments/u0.01875"

python2 do_power_analysis.py
python2 create_AFS_figures.py
python2 create_number_of_divisions_table.py > S1table.txt
python2 create_Splots_figures.py ../model/experiments/u0.01875/1_0_0_outs_10/Mar1pipe_out_Fri_Mar__2_20_42_31_2018/ ../model/experiments/u0.01875/1_0_02_outs_10/Mar1pipe_out_Fri_Mar__2_21_18_46_2018/ "10"
python2 create_cluster_advantage_plots.py
python2 create_fanning_plots.py
python2 create_mixing_analysis.py "../model/experiments/u0.01875/1_0_0" "10" "Noturnover"
python2 create_mixing_analysis.py "../model/experiments/u0.01875/1_0_005" "10" "Turnover005"
python2 create_mixing_analysis.py "../model/experiments/u0.01875/1_0_01" "10" "Turnover01"
python2 create_mixing_analysis.py "../model/experiments/u0.01875/1_0_02" "10" "Turnover02"
python2 create_mixing_analysis.py "../model/experiments/u0.01875/1_0_065" "10" "Turnover065"
python2 create_trees.py "../model/experiments/u0.01/0_0_065" "10" "Turnover065"
python2 create_trees.py "../model/experiments/u0.01/0_0_02" "10" "Turnover02"
python2 create_trees.py "../model/experiments/u0.01/0_0_01" "10" "Turnover01"
python2 create_trees.py "../model/experiments/u0.01/0_0_005" "10" "Turnover005"
python2 create_trees.py "../model/experiments/u0.01/0_0_0" "10" "NoTurnover"


