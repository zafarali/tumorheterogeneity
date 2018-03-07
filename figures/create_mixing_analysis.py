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
ALL_SEEDS = ['10','100','102','15','3','3318','33181','33185','33186','34201810','342018101','342018102','8','9','99']
import json
from analysis.newick import loads
from scipy.stats import spearmanr
from decimal import Decimal

"""
Ranked Distance plots
"""


def plot_ranked_mixing_analysis(root_folder, label, fig, to_plot, color='r',
                                sort=True, total_rows=1, row_to_start=0, show_labels=False):
    search_folder = root_folder + '*/Mar*/mixing_analysis.npy'
    print('Search foler:', search_folder)
    mixing_files = glob(search_folder)

    all_data = []
    print(mixing_files)

    for file_ in mixing_files:
        data = np.load(file_)
        data = data[np.argsort(data[:, 0, 0])]  # sort the data according to how far from center.

        if fig is None:
            fig = plt.figure(figsize=(13, 3))
        max_value = 0
        all_axes = []
        for idx, i in enumerate([0, 1, 2, 3, -4, -3, -2, -1]):  # TODO: for now just enumerate through the extremes
            ax = fig.add_subplot(total_rows, 8, (8 * row_to_start) + idx + 1)

            distance_to_tumor_COM, distance_to_cluster_COM, p_sgas, n_sgas, p_pgas, n_pgas, cum_S = data[i]
            distance_to_tumor_COM = distance_to_tumor_COM[0]

            p_pgas = np.array(p_pgas)[1:]
            p_sgas = np.array(p_sgas)[1:]
            n_pgas = np.array(n_pgas)[1:]
            n_sgas = np.array(n_sgas)[1:]
            cum_S = np.array(cum_S)[1:]
            distance_to_cluster_COM = np.array(distance_to_cluster_COM)[1:]

            if sort:
                sorted_distances = np.argsort(distance_to_cluster_COM)
                p_pgas = p_pgas[sorted_distances]
                n_pgas = n_pgas[sorted_distances]
                n_sgas = n_sgas[sorted_distances]
                p_sgas = p_sgas[sorted_distances]
                cum_S = cum_S[sorted_distances]
                distances = np.arange(0, n_sgas.shape[0], 1)
            else:
                distances = distance_to_cluster_COM

            if to_plot == 'n_pgas':
                ax.scatter(distances, n_pgas, \
                           label='{}'.format(label), color=color)
                rho, p = spearmanr(distances, n_pgas)
                max_value = np.max([np.max(n_pgas), max_value])

            elif to_plot == 'p_pgas':
                ax.scatter(distances, p_pgas, \
                           label='{}'.format(label, distance_to_tumor_COM), color=color)
                rho, p = spearmanr(distances, p_pgas)
                max_value = np.max([np.max(p_pgas), max_value])


            elif to_plot == 'n_sgas':
                ax.scatter(distances, n_sgas, \
                           label='{}'.format(label, distance_to_tumor_COM), color=color)
                rho, p = spearmanr(distances, n_sgas)
                max_value = np.max([np.max(n_sgas), max_value])


            elif to_plot == 'p_sgas':
                ax.scatter(distances, p_sgas, \
                           label='{}'.format(label, distance_to_tumor_COM), color=color)
                rho, p = spearmanr(distances, p_sgas)
                max_value = np.max([np.max(p_sgas), max_value])

            elif to_plot == 'cum_S':
                ax.plot(distances, cum_S, \
                        label='{}'.format(label, distance_to_tumor_COM), color=color)
                max_value = np.max([np.max(cum_S), max_value])
                rho, p = spearmanr(distances, cum_S)
            if idx == 0: ax.set_ylabel(label+'\n'+to_plot)
            if show_labels: ax.set_xlabel('Distance from \ncluster COM')
            ax.set_title('x={:3g},\nrho={:.2g},\np={:.2e}'.format(distance_to_tumor_COM, rho, Decimal(p)))
            all_axes.append(ax)
            sns.despine()
        # ax.legend()
        for ax in all_axes:
            ax.set_ylim([0, max])

        break


PATH, SEED = sys.argv[1], sys.argv[2]
thing_to_plot = 'n_pgas'
fig = plt.figure(figsize=(14, 10))
plot_ranked_mixing_analysis(PATH+'0_outs_'+SEED, 'No Turnover', fig, thing_to_plot, sort=False, total_rows=5, row_to_start=0)
plot_ranked_mixing_analysis(PATH+'005_outs_'+SEED, 'Turnover (d=0.05)', fig, thing_to_plot, sort=False, total_rows=5, row_to_start=1)
plot_ranked_mixing_analysis(PATH+'01_outs_'+SEED, 'Turnover (d=0.1)', fig, thing_to_plot, sort=False, total_rows=5, row_to_start=2)
plot_ranked_mixing_analysis(PATH+'02_outs_'+SEED, 'Turnover (d=0.2)', fig, thing_to_plot, sort=False, total_rows=5, row_to_start=3)
plot_ranked_mixing_analysis(PATH+'065_outs_'+SEED, 'Turnover (d=0.65)', fig, thing_to_plot, sort=False, total_rows=5, row_to_start=4, show_labels=True)
fig.tight_layout(h_pad=1)

