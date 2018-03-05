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
python2 do_big_power_analysis.py
python2 create_AFS_figures.py
python2 create_number_of_divisions_table.py > S1table.txt
python2 create_SvsN_figures.py
python2 create_Splots_figures.py
python2 create_cluster_advantage_plots.py
python2 create power_analysis_big.py
python2 create_fanning_plots.py
python2 create_mixing_analysis.py
python2 create_mixing_analysis.py "../model/experiments/u0.01875/1_0_0" "10" "No turnover"
python2 create_mixing_analysis.py "../model/experiments/u0.01875/1_0_005_outs" "10" "Turnover005"
python2 create_mixing_analysis.py "../model/experiments/u0.01875/1_0_01_outs" "10" "Turnover01"
python2 create_mixing_analysis.py "../model/experiments/u0.01875/1_0_02_outs" "10" "Turnover02"
python2 create_mixing_analysis.py "../model/experiments/u0.01875/1_0_065_outs" "10" "Turnover065"
python2 create_trees.py "../model/experiments/u0.01875/1_0_065" "10" > turnover065_trees_10.txt
python2 create_trees.py "../model/experiments/u0.01875/1_0_02" "10" > turnover02_trees_10.txt
python2 create_trees.py "../model/experiments/u0.01875/1_0_01" "10" > turnover01_trees_10.txt
python2 create_trees.py "../model/experiments/u0.01875/1_0_005" "10" > turnover005_trees_10.txt
python2 create_trees.py "../model/experiments/u0.01875/1_0_0" "10" > no_turnover_trees_10.txt


