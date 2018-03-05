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
ALL_SEEDS = ['10', '102', '100', '15', '3', '3318', '33181', '33185', '9', '99', '33186']
import json
from analysis.newick import loads

"""
Fanning Plots
"""

# plots the "fanning" plots
def plot_diffs(root_folder, seeds, k=0.01, pts=100, cutoff='00'):
    fig = plt.figure(figsize=(13, 3))

    ax = fig.add_subplot(131)
    folder = root_folder + '/1_0'
    d2 = '0'
    No_turnover = data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d=d2)
    Turnover = data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d='005')
    to_plot = {}
    for k in No_turnover.keys():
        avg = np.mean(k)
        to_plot[k] = [Turnover[(1,)][0], avg * (Turnover[(1,)][1] - No_turnover[(1,)][1]) + No_turnover[k][1]]

    plot_it(to_plot, ax)
    plot_it(Turnover, ax, '--')
    ax.set_title('d=0.05')
    ax.set_xlabel('Distance from COM of Tumor')
    ax.set_ylabel('S')
    Turnover = data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d='01')

    ax = fig.add_subplot(132)
    to_plot = {}
    for k in No_turnover.keys():
        avg = np.mean(k)
        to_plot[k] = [Turnover[(1,)][0], avg * (Turnover[(1,)][1] - No_turnover[(1,)][1]) + No_turnover[k][1]]

    plot_it(to_plot, ax)
    plot_it(Turnover, ax, '--')
    ax.set_xlabel('Distance from COM of Tumor')
    ax.set_title('d=0.1')
    ax.set_ylabel('S')
    ax = fig.add_subplot(133)
    Turnover = data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d='02')

    to_plot = {}
    for k in No_turnover.keys():
        avg = np.mean(k)
        to_plot[k] = [Turnover[(1,)][0], avg * (Turnover[(1,)][1] - No_turnover[(1,)][1]) + No_turnover[k][1]]

    plot_it(to_plot, ax)
    plot_it(Turnover, ax, '--')
    ax.set_title('d=0.2')
    ax.set_xlabel('Distance from COM of Tumor')
    ax.set_ylabel('S')
    fig.tight_layout(h_pad=1)
    return fig

plot_diffs('../model/experiments/u0.01875/',ALL_SEEDS).savefig('./Splot-fanning.pdf')

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

plot_ranked_mixing_analysis('../model/experiments/u0.01875/1_0_0_outs_10', 'No Turnover', None, 'n_pgas')
plot_ranked_mixing_analysis('../model/experiments/u0.01875/1_0_0_outs_10', 'No Turnover', None, 'cum_S')
plot_ranked_mixing_analysis('../model/experiments/u0.01875/1_0_065_outs_10', 'Turnover 065', None, 'n_pgas')
plot_ranked_mixing_analysis('../model/experiments/u0.01875/1_0_065_outs_10', 'Turnover 065', None, 'cum_S')
