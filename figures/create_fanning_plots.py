import sys
sys.path.append('..')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from glob import glob
sns.set_style('white')
sns.set_context('paper', font_scale=1.5)
from figtools import *
ALL_SEEDS = ['10','100','102','15','3','3318','33181','33185','33186','34201810','342018101','342018102','8','9','99']

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



def plot_relative_increases(root_folder, seeds):
    fig = plt.figure(figsize=(8, 3))

    ax = fig.add_subplot(111)
    folder = root_folder + '/1_0'

    No_turnover = data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d='0', max_dist=300)
    Turnover005 = data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d='005', max_dist=300)
    Turnover01 = data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d='01', max_dist=300)
    Turnover02 = data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d='02', max_dist=300)
    Turnover065 = data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d='065', max_dist=300)

    
    x,y = [No_turnover[(1,)][0], (No_turnover[(1,)][1] - No_turnover[(1,)][1])/(No_turnover[(1,)][1]) ]

    ax.plot(x, y, label='d=0')
    x,y = [No_turnover[(1,)][0], (Turnover005[(1,)][1] - No_turnover[(1,)][1])/(No_turnover[(1,)][1]) ]

    ax.plot(x, y, label='d=0.05')
    to_plot['01'] = [No_turnover[(1,)][0], (Turnover01[(1,)][1] - No_turnover[(1,)][1])/(No_turnover[(1,)][1]) ]

    ax.plot(x, y, label='d=0.1')
    to_plot['02'] = [No_turnover[(1,)][0], (Turnover02[(1,)][1] - No_turnover[(1,)][1])/(No_turnover[(1,)][1]) ]

    ax.plot(x, y, label='d=0.2')
    to_plot['065'] = [No_turnover[(1,)][0], (Turnover065[(1,)][1] - No_turnover[(1,)][1])/(No_turnover[(1,)][1]) ]

    ax.plot(x, y, label='d=0.65')

    ax.set_xlabel('Distance from COM of Tumor')
    ax.set_ylabel('Relative Increase in S(1)')
    ax.legend()
    fig.tight_layout(h_pad=1)
    return fig


plot_relative_increases('../model/experiments/u0.01875/',ALL_SEEDS).savefig('./Splot-relative-increase.pdf')

