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
N = 10**8
b = 0.69
R_f = lambda d: (3 * N / (4. * np.pi * (1 - d / b))) ** (1.0 / 3.0)

# plots the "fanning" plots
def plot_diffs(root_folder, seeds, k=0.01, pts=100, cutoff='00', gamma=0.5, mu=0.02):
    fig = plt.figure(figsize=(10, 3))

    ax = fig.add_subplot(131)
    folder = root_folder + '/0_0'
    d2 = '0'
    No_turnover = data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d=d2)
    Turnover = data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d='005')
    to_plot = {}
    d = 0.05
    for k in No_turnover.keys():
        avg = np.mean(k)
        to_plot[k] = [Turnover[(1,)][0], avg*mu*d*(R_f(d) - Turnover[(1,)][0])/(gamma * b) + (1+d) * No_turnover[k][1]]

    plot_it(to_plot, ax, '--')
    plot_it(Turnover, ax)
    ax.plot(No_turnover[(1,)][0], No_turnover[(1,)][1], ':k')
    ax.set_title('d=0.05')
    ax.set_xlabel('Distance from COM of Tumor')
    ax.set_ylabel('S')
    Turnover = data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d='01')

    ax = fig.add_subplot(132)
    to_plot = {}
    d = 0.1
    for k in No_turnover.keys():
        avg = np.mean(k)
        to_plot[k] = [Turnover[(1,)][0], avg*mu*d*(R_f(d) - Turnover[(1,)][0])/(gamma * b) + (1+d) * No_turnover[k][1]]

    plot_it(to_plot, ax, '--')
    plot_it(Turnover, ax)
    ax.plot(No_turnover[(1,)][0], No_turnover[(1,)][1], ':k')
    ax.set_xlabel('Distance from COM of Tumor')
    ax.set_title('d=0.1')
    ax.set_ylabel('S')
    ax = fig.add_subplot(133)
    Turnover = data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d='02')

    to_plot = {}
    d=0.2
    for k in No_turnover.keys():
        avg = np.mean(k)
        to_plot[k] = [Turnover[(1,)][0], avg*mu*d*(R_f(d) - Turnover[(1,)][0])/(gamma * b) + (1+d) * No_turnover[k][1]]

    plot_it(to_plot, ax, '--')
    plot_it(Turnover, ax)
    ax.plot(No_turnover[(1,)][0], No_turnover[(1,)][1], ':k')
    ax.set_title('d=0.2')
    ax.set_xlabel('Distance from COM of Tumor')
    ax.set_ylabel('S')
    from matplotlib.lines import Line2D

    custom_lines = [Line2D([0], [0], color='k', linestyle='--'),
                    Line2D([0], [0], color='k', linestyle='-'),
                    Line2D([0], [0], color='k', linestyle=':')]
    ax.legend(custom_lines, ['Prediction', 'Observed', 'S(1), d=0'])
    fig.tight_layout(h_pad=1)
    return fig

plot_diffs('../model/experiments/u0.01/',ALL_SEEDS, mu=0.01).savefig('./Splot-fanning-001.pdf')
# plot_diffs('../model/experiments/u0.01/',ALL_SEEDS).savefig('./Splot-fanning-002.pdf')

