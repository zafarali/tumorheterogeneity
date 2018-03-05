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
def plot_ranked_mixing_analysis(root_folder, label, ax, to_plot, color='r'):
    mixing_files = glob(root_folder+'*/Mar*/mixing_analysis.npy')

    all_data = []
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

            sorted_distances = np.argsort(distance_to_cluster_COM)
            ranked_p_pgas = p_pgas[sorted_distances]
            ranked_n_pgas = n_pgas[sorted_distances]
            ranked_n_sgas = n_sgas[sorted_distances]
            ranked_p_sgas = p_sgas[sorted_distances]
            ranked_cum_S = cum_S[sorted_distances]

            if to_plot == 'n_pgas':
                ax.scatter(np.arange(0, ranked_n_pgas.shape[0], 1), ranked_n_pgas, \
                            label='{}'.format(label), color=color)
                ax.set_ylabel('# Private GAS')
                ax.set_title(' distance: {:3g}'.format(distance_to_tumor_COM))
            elif to_plot == 'p_pgas':
                ax.scatter(np.arange(0, ranked_p_pgas.shape[0], 1), ranked_p_pgas, \
                            label='{}'.format(label, distance_to_tumor_COM), color=color)
                ax.set_ylabel('% Private GAS')
                ax.set_title(' distance: {:3g}'.format(distance_to_tumor_COM))

            elif to_plot == 'n_sgas':
                ax.scatter(np.arange(0, ranked_n_sgas.shape[0], 1), ranked_n_sgas, \
                            label='{}'.format(label, distance_to_tumor_COM), color=color)
                ax.set_ylabel('# Shared GAS')
                ax.set_title(' distance: {:3g}'.format(distance_to_tumor_COM))

            elif to_plot == 'p_sgas':
                ax.scatter(np.arange(0, ranked_p_sgas.shape[0], 1), ranked_p_sgas, \
                            label='{}'.format(label, distance_to_tumor_COM), color=color)
                ax.set_ylabel('% Shared GAS')
                ax.set_title(' distance: {:3g}'.format(distance_to_tumor_COM))

            elif to_plot == 'cum_S':
                ax.plot(np.arange(0, ranked_cum_S.shape[0], 1), ranked_cum_S, \
                            label='{}'.format(label, distance_to_tumor_COM), color=color)
                ax.set_ylabel('S(n)')
                ax.set_title(' distance: {:3g}'.format(distance_to_tumor_COM))


            ax.set_xlabel('Ranked Distance \nfrom Cluster center')
            ax.legend()
        fig.tight_layout(h_pad=1)
        fig.savefig('./ranked_mixing_{}_{}.pdf'.format(label, to_plot))
        break


PATH, SEED, title = sys.argv[1], sys.argv[2], sys.argv[3]

plot_ranked_mixing_analysis(PATH+'_outs_'+SEED, title, None, 'n_pgas')
plot_ranked_mixing_analysis(PATH+'_outs_'+SEED, title, None, 'p_pgas')
plot_ranked_mixing_analysis(PATH+'_outs_'+SEED, title, None, 'cum_S')