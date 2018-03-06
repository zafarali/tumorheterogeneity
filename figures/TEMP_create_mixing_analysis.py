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


"""
Ranked Distance plots
"""
def plot_ranked_mixing_analysis(root_folder, label, ax, to_plot, color='r', sort=True):
    search_folder = root_folder+'*/MarMixingAnalysisOnly*/mixing_analysis.npy'
    print('Search foler:', search_folder)
    mixing_files = glob(search_folder)

    all_data = []
    print(mixing_files)

    for file_ in mixing_files:
        data = np.load(file_)
        data = data[np.argsort(data[:, 0, 0])] # sort the data according to how far from center.

        fig = plt.figure(figsize=(15, 8))
        subplots = '24'

        for idx, i in enumerate([0, 1, 2, 3, -4, -3, -2, -1]): #TODO: for now just enumerate through the extremes
            ax = fig.add_subplot(subplots+str(idx+1))
            distance_to_tumor_COM, distance_to_cluster_COM, p_sgas, n_sgas, p_pgas, n_pgas, cum_S = data[i]
            distance_to_tumor_COM = distance_to_tumor_COM[0]

            p_pgas = np.array(p_pgas)
            p_sgas = np.array(p_sgas)
            n_pgas = np.array(n_pgas)
            n_sgas = np.array(n_sgas)
            cum_S = np.array(cum_S)

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
                ax.set_ylabel('# Private GAS')
                ax.set_title(' distance: {:3g}'.format(distance_to_tumor_COM))
            elif to_plot == 'p_pgas':
                ax.scatter(distances, p_pgas, \
                            label='{}'.format(label, distance_to_tumor_COM), color=color)
                ax.set_ylabel('% Private GAS')
                ax.set_title(' distance: {:3g}'.format(distance_to_tumor_COM))

            elif to_plot == 'n_sgas':
                ax.scatter(distances, n_sgas, \
                            label='{}'.format(label, distance_to_tumor_COM), color=color)
                ax.set_ylabel('# Shared GAS')
                ax.set_title(' distance: {:3g}'.format(distance_to_tumor_COM))

            elif to_plot == 'p_sgas':
                ax.scatter(distances, p_sgas, \
                            label='{}'.format(label, distance_to_tumor_COM), color=color)
                ax.set_ylabel('% Shared GAS')
                ax.set_title(' distance: {:3g}'.format(distance_to_tumor_COM))

            elif to_plot == 'cum_S':
                ax.plot(distances, cum_S, \
                            label='{}'.format(label, distance_to_tumor_COM), color=color)
                ax.set_ylabel('S(n)')
                ax.set_title(' distance: {:3g}'.format(distance_to_tumor_COM))

            if sort:
                ax.set_xlabel('Ranked Distance \nfrom Cluster center')
            else:
                ax.set_xlabel('Distance from \nCluster Center')
            ax.legend()
        fig.tight_layout(h_pad=1)
        fig.savefig(root_folder+'./ranked_mixing_{}_{}_{}.pdf'.format(label, to_plot, sort))
        break


PATH, title = sys.argv[1], sys.argv[2]

plot_ranked_mixing_analysis(PATH, title, None, 'n_pgas', sort=True)
plot_ranked_mixing_analysis(PATH, title, None, 'p_pgas', sort=True)
plot_ranked_mixing_analysis(PATH, title, None, 'cum_S', sort=True)
plot_ranked_mixing_analysis(PATH, title, None, 'n_pgas', sort=False)
plot_ranked_mixing_analysis(PATH, title, None, 'p_pgas', sort=False)
plot_ranked_mixing_analysis(PATH, title, None, 'cum_S', sort=False)