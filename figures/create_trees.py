import sys
sys.path.append('..')
import random
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from mpl_toolkits.axes_grid.inset_locator import inset_axes
sns.set_style('white')
sns.set_context('paper', font_scale=1.5)
from glob import glob
from figtools import *
from analysis.newick_tree_maker import neighbor_joining_tree
import json
from analysis.newick import loads

"""
Trees
"""

def create_trees(root_folder, seed):
    cluster_search = root_folder +'_outs_'+seed+'/Mar*/cluster_data.json'
    mixing_search = root_folder +'_outs_'+seed+'/Mar*/mixing_analysis.npy'
    print('looking in:',cluster_search, mixing_search)
    cluster_data = glob(cluster_search)
    mixing_data = glob(mixing_search)
    for cluster_data_,mixing_data_ in zip(cluster_data, mixing_data):
        data = json.load(open(cluster_data_, 'r'))
        mixing_data = np.load(mixing_data_)
        sorted_idx = np.argsort(mixing_data[:, 0, 0])
        for idx in sorted_idx[:4].tolist() + sorted_idx[-4:].tolist():
            print(idx)
            print(data[idx])
            distance_from_COM = np.sqrt(np.sum(np.array(data[idx]['COM'])**2))
            tree_string, mapper = neighbor_joining_tree(data[idx]['genotypes'])
            print('START'*20)
            print(mapper)
            print(distance_from_COM)
            print(tree_string)
            print(loads(tree_string)[0].ascii_art().encode('utf-8'))
            print('END'*40)

PATH, SEED = sys.argv[1], sys.argv[2]
create_trees(PATH, SEED)
# create_trees('../model/no_death', '1')


