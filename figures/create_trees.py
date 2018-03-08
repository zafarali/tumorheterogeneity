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
def create_trees(root_folder, seed, title, fig, total_rows=1, row_to_start=0):
    cluster_search = root_folder +'_outs_'+seed+'/Mar*/cluster_data.json'
    mixing_search = root_folder +'_outs_'+seed+'/Mar*/mixing_analysis.npy'
    print('looking in:',cluster_search, mixing_search)
    cluster_data = glob(cluster_search)
    mixing_data = glob(mixing_search)



    for z, (cluster_data_,mixing_data_) in enumerate(zip(cluster_data, mixing_data)):
        data = json.load(open(cluster_data_, 'r'))
        # mixing_data_ = np.load(mixing_data_)
        distances_ = [ np.sqrt(np.sum(np.array(data_details['COM'])**2)) for data_details in data ]
        sorted_idx = np.argsort(distances_) # get the distances ranked

        if fig is None:
            fig = plt.figure(figsize=(13, 3))
        all_axes = []
        for panel_count, idx in enumerate(sorted_idx[:4].tolist() + sorted_idx[-4:].tolist()):
            ax = fig.add_subplot(total_rows, 8, (8 * row_to_start) + panel_count + 1)
            if panel_count == 0: ax.set_ylabel(title)
            distance_from_COM = np.sqrt(np.sum(np.array(data[idx]['COM'])**2))
            data_, all_snps = prepare_snp_data(data[idx]['genotypes'])

            tree = Node(data_, all_snps=all_snps.most_common())
            tree.resolve()
            tree.plot(1, 0, ax=ax)
            ax.set_title('x={:3g}'.format(distance_from_COM))
            sns.despine(ax=ax, top=True, bottom=True, left=False, right=True)
            ax.tick_params(
                axis='x',  # changes apply to the x-axis
                which='both',  # both major and minor ticks are affected
                bottom='off',  # ticks along the bottom edge are off
                top='off',  # ticks along the top edge are off
                labelbottom='off')
            all_axes.append(ax)

        max_ymax = 0
        min_ymin = -100
        for ax in all_axes:
            ymin, ymax = ax.get_ylim()
            min_ymin = np.min([ymin, min_ymin])
            max_ymax = np.max([ymax, max_ymax])
        for ax in all_axes:
            ax.set_ylim([min_ymin, max_ymax])

    return fig

PATH, SEED = sys.argv[1], sys.argv[2]
# create_trees('../model/no_death', '1')
fig = plt.figure(figsize=(14,10))
create_trees(PATH+'/1_0_0', SEED, 'No Turnover', fig, total_rows=5, row_to_start=0)
create_trees(PATH+'/1_0_005', SEED, 'Turnover (d=0.05)', fig, total_rows=5, row_to_start=1)
create_trees(PATH+'/1_0_01', SEED, 'Turnover (d=0.1)', fig, total_rows=5, row_to_start=2)
create_trees(PATH+'/1_0_02', SEED, 'Turnover (d=0.2)', fig, total_rows=5, row_to_start=3)
create_trees(PATH+'/1_0_065', SEED, 'Turnover (d=0.65)', fig, total_rows=5, row_to_start=4)
fig.tight_layout(h_pad=1)
fig.savefig('./trees_selection.pdf')

fig = plt.figure(figsize=(14,10))
create_trees(PATH+'/0_0_0', SEED, 'No Turnover', fig, total_rows=5, row_to_start=0)
create_trees(PATH+'/0_0_005', SEED, 'Turnover (d=0.05)', fig, total_rows=5, row_to_start=1)
create_trees(PATH+'/0_0_01', SEED, 'Turnover (d=0.1)', fig, total_rows=5, row_to_start=2)
create_trees(PATH+'/0_0_02', SEED, 'Turnover (d=0.2)', fig, total_rows=5, row_to_start=3)
create_trees(PATH+'/0_0_065', SEED, 'Turnover (d=0.65)', fig, total_rows=5, row_to_start=4)
fig.tight_layout(h_pad=1)
fig.savefig('./trees_no_selection.pdf')

fig = plt.figure(figsize=(14,10))
create_trees(PATH+'/10_0_0', SEED, 'No Turnover', fig, total_rows=5, row_to_start=0)
create_trees(PATH+'/10_0_005', SEED, 'Turnover (d=0.05)', fig, total_rows=5, row_to_start=1)
create_trees(PATH+'/10_0_01', SEED, 'Turnover (d=0.1)', fig, total_rows=5, row_to_start=2)
create_trees(PATH+'/10_0_02', SEED, 'Turnover (d=0.2)', fig, total_rows=5, row_to_start=3)
create_trees(PATH+'/10_0_065', SEED, 'Turnover (d=0.65)', fig, total_rows=5, row_to_start=4)
fig.tight_layout(h_pad=1)
fig.savefig('./trees_extreme_selection.pdf')

