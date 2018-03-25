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




"""
No selection case.
"""
def plot_relative_increases(root_folder, seeds, d_append=''):
    fig = plt.figure(figsize=(8, 3))

    ax = fig.add_subplot(111)
    folder = root_folder + '/0_0'

    No_turnover = data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d='0'+d_append, max_dist=300)
    Turnover005 = data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d='005'+d_append, max_dist=300)
    Turnover01 = data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d='01'+d_append, max_dist=300)
    Turnover02 = data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d='02'+d_append, max_dist=300)
    Turnover065 = data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d='065'+d_append, max_dist=300)


    x,y = [No_turnover[(1,)][0], (Turnover005[(1,)][1] - No_turnover[(1,)][1])/(No_turnover[(1,)][1]) ]

    ax.plot(x, y, label='d=0.05')
    x,y = [No_turnover[(1,)][0], (Turnover01[(1,)][1] - No_turnover[(1,)][1])/(No_turnover[(1,)][1]) ]

    ax.plot(x, y, label='d=0.1')
    x,y = [No_turnover[(1,)][0], (Turnover02[(1,)][1] - No_turnover[(1,)][1])/(No_turnover[(1,)][1]) ]

    ax.plot(x, y, label='d=0.2')
    x,y = [No_turnover[(1,)][0], (Turnover065[(1,)][1] - No_turnover[(1,)][1])/(No_turnover[(1,)][1]) ]

    ax.plot(x, y, label='d=0.65')

    ax.set_xlabel('Distance from COM of Tumor')
    ax.set_ylabel('Relative Increase in S(1)')
    ax.legend()
    ax.semilogy(basey=10)
    fig.tight_layout(h_pad=1)
    return fig


def S_plot_no_selection(root_folder, seeds, d_append=''):
    fig = plt.figure(figsize=(13, 3))


    d = '005'
    ax = fig.add_subplot(141)

    cmap_ax = ax  # for later use

    folder = root_folder + '/0_0'

    plot_it(data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d=d+d_append), ax)
    #     plot_it(data_to_plot(folder, seeds, yaxis='S_list'+cutoff+'_big', mode=2, d=d), ax)
    plot_density(ax, folder, seeds, d=d+d_append)

    ax.set_title('(b) Surface Turnover, d=0.' + d[1:])
    ax.set_xlim([50, 325])
    ax.set_ylabel('Mean S(n)')
    ax.set_xlabel('Distance from Centre \n of Tumor (cells)')
    #     ax.set_ylim([0,1])

    d = '01'
    ax = fig.add_subplot(142)
    folder = root_folder + '/0_0'
    plot_it(data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d=d+d_append), ax)
    #     plot_it(data_to_plot(folder, seeds, yaxis='S_list'+cutoff+'_big', mode=2, d=d), ax)
    plot_density(ax, folder, seeds, d=d+d_append)

    ax.set_title('(c) Surface Turnover, d=0.' + d[1:])
    ax.set_ylabel('Mean S(n)')
    x_max = 350 if d == '02' else 325
    x_max = 750 if d == '065' else x_max
    ax.set_xlim([50, x_max])
    ax.set_xlabel('Distance from Centre \n of Tumor (cells)')

    d = '02'
    ax = fig.add_subplot(143)
    folder = root_folder + '/0_0'

    plot_it(data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d=d+d_append), ax)
    #     plot_it(data_to_plot(folder, seeds, yaxis='S_list'+cutoff+'_big', mode=2, d=d), ax)
    plot_density(ax, folder, seeds, d=d+d_append)

    ax.set_title('(d) Surface Turnover, d=0.' + d[1:])
    ax.set_ylabel('Mean S(n)')
    x_max = 350 if d == '02' else 325
    x_max = 750 if d == '065' else x_max
    ax.set_xlim([50, x_max])
    ax.set_xlabel('Distance from Centre \n of Tumor (cells)')
    #     ax.set_ylim([0,50])

    d = '065'
    ax = fig.add_subplot(144)
    folder = root_folder + '/0_0'

    plot_it(data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d=d+d_append), ax)
    #     plot_it(data_to_plot(folder, seeds, yaxis='S_list'+cutoff+'_big', mode=2, d=d), ax)
    plot_density(ax, folder, seeds, d=d+d_append)

    ax.set_title('(e) Surface Turnover, d=0.' + d[1:])
    ax.set_ylabel('Mean S(n)')
    x_max = 350 if d == '02' else 325
    #     x_max = 750 if d=='065' else x_max
    ax.set_xlim([50, x_max])
    ax.set_xlabel('Distance from Centre \n of Tumor (cells)')
    #     ax.set_ylim([0,50])

    fig.tight_layout(pad=1, w_pad=0.5)

    subaxes = inset_axes(cmap_ax,
                         width="90%",  # width = 30% of parent_bbox
                         height=0.15,  # height : 1 inch
                         loc=2)

    create_colorbar(subaxes, with_yellow=True)
    return fig

plot_relative_increases('../model/experiments/u0.01/',ALL_SEEDS).savefig('./Splot-relative-increase.pdf')

SMALL_SEEDS = ['1', '2']

fig = plot_relative_increases('../model/experiments/u0.01small/',SMALL_SEEDS, d_append='_107')
fig.suptitle('10^7')
fig.savefig('./Splot-relative-increase107.pdf')

fig = S_plot_no_selection('../model/experiments/u0.01small/',SMALL_SEEDS, d_append='_107')
fig.suptitle('10^7')
fig.savefig('./Splot-107.pdf')


fig = plot_relative_increases('../model/experiments/u0.01small/',SMALL_SEEDS, d_append='_106')
fig.suptitle('10^6')
fig.savefig('./Splot-relative-increase106.pdf')

fig = S_plot_no_selection('../model/experiments/u0.01small/',SMALL_SEEDS, d_append='_106')
fig.suptitle('10^6')
fig.savefig('./Splot-106.pdf')


fig = plot_relative_increases('../model/experiments/u0.01small/',SMALL_SEEDS, d_append='_105')
fig.suptitle('10^5')
fig.savefig('./Splot-relative-increase105.pdf')


fig = S_plot_no_selection('../model/experiments/u0.01small/',SMALL_SEEDS, d_append='_105')
fig.suptitle('10^5')
fig.savefig('./Splot-105.pdf')


fig = plot_relative_increases('../model/experiments/u0.01small/',SMALL_SEEDS, d_append='_104')
fig.suptitle('10^4')
fig.savefig('./Splot-relative-increase104.pdf')

fig = S_plot_no_selection('../model/experiments/u0.01small/',SMALL_SEEDS, d_append='_104')
fig.suptitle('10^4')
fig.savefig('./Splot-104.pdf')