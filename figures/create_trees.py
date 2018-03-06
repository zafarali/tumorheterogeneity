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
from analysis.newick_tree_maker import prepare_snp_data, Node
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
    for z, (cluster_data_,mixing_data_) in enumerate(zip(cluster_data, mixing_data)):
        data = json.load(open(cluster_data_, 'r'))
        mixing_data_ = np.load(mixing_data_)
        sorted_idx = np.argsort(mixing_data_[:, 0, 0])
        panel = '24'
        fig = plt.figure(figsize=(18, 10))

        for panel_count, idx in enumerate(sorted_idx[:4].tolist() + sorted_idx[-4:].tolist()):
            panel_ = panel + str(panel_count+1)
            ax = fig.add_subplot(panel_)
            # print(data[idx])
            distance_from_COM = np.sqrt(np.sum(np.array(data[idx]['COM'])**2))
            data_, all_snps = prepare_snp_data(data[idx]['genotypes'])

            tree = Node(data_, all_snps=all_snps)
            tree.resolve()
            tree.plot(1, 0, ax=ax)
            ax.set_title('distance:{:3g}'.format(distance_from_COM))

        fig.savefig(root_folder+'_outs_'+seed+'/tree_plots_{}.pdf'.format(z))

# PATH, SEED = sys.argv[1], sys.argv[2]
create_trees(PATH, SEED)
# create_trees('../model/no_death', '1')


