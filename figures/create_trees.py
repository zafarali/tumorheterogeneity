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

def create_trees(root_folder, seeds):
    for seed in seeds:
        folder_search = root_folder +'_outs_'+seed+'/Mar*/cluster_data.json'
        print(folder_search)
        cluster_data = glob(folder_search)
        print(cluster_data)
        for cluster_data_ in cluster_data:
            data = json.load(open(cluster_data_, 'r'))
            for data_ in data:
                distance_from_COM = np.sqrt(np.sum(np.array(data_['COM'])**2))
                tree_string, mapper = neighbor_joining_tree(data_['genotypes'])
                print('START'*20)
                print(mapper)
                print(distance_from_COM)
                print(tree_string)
                print(loads(tree_string)[0].ascii_art().encode('utf-8'))
                print('END'*40)

PATH, SEED = sys.argv[1], sys.argv[1]
create_trees(PATH, SEED)


