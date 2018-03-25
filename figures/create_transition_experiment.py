import sys
import random
import statsmodels.api as sm
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
import json
from figtools import *
from decimal import Decimal

sns.set_style('white')
sns.set_context('paper', font_scale=1.5)
from glob import glob
YELLOW = sns.xkcd_rgb["amber"]
params = {'mathtext.default': 'regular' }
plt.rcParams.update(params)
CTC_SIZES = [(2,7), (8,12), (13,17), (18,22),(23,30)]

ALTERNATE_COLORS = sns.xkcd_palette(["vomit green", "orangish brown", "azure", "dark pink"])
small_c, med_c, big_c, biggest_c = ALTERNATE_COLORS


def S_plot_transition(root_folder, seeds):
    fig = plt.figure(figsize=(13, 6))


    """
    Small Tumor of Size 10^6
    """
    #
    # d_append = '_106'
    # folder = root_folder + '/0_0'
    # for i, d in enumerate(['0', '005', '01', '02', '065']):
    #     ax = fig.add_subplot(2,5,i+1)
    #     plot_it(data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d=d+d_append), ax)
    #     plot_density(ax, folder, seeds, d=d+d_append)
    #     ax.set_title('$N=10^6, d=0.' + d[1:]+'$')
    #     ax.set_xlim([20, 325])
    #     ax.set_ylim([0, 30])
    #     ax.set_ylabel('Mean S(n)')
    #     ax.set_xlabel('Distance from Centre \n of Tumor (cells)')
    #
    #     if i == 0: cmap_ax = ax  # for later use
    #
    #
    # """
    # Large Tumor of Size 10^8 but no new mutations are introduced after 10^6.
    # """
    # d_append = '_108'
    # for i, d in enumerate(['0', '005', '01', '02', '065']):
    #     ax = fig.add_subplot(2,5,6+i)
    #     plot_it(data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d=d+d_append), ax)
    #     plot_density(ax, folder, seeds, d=d+d_append)
    #     ax.set_title('$N=10^8, d=0.' + d[1:]+'$')
    #     ax.set_xlim([20, 325])
    #     ax.set_ylim([0, 30])
    #     ax.set_ylabel('Mean S(n)')

    d='01'
    d_append = '_106'
    folder = root_folder + '/0_0'

    ax = fig.add_subplot(1, 2, 1)
    plot_it(data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d=d + d_append), ax)
    plot_density(ax, folder, seeds, d=d + d_append)
    ax.set_title('(a) $N=10^6, d=0.' + d[1:] + '$')
    ax.set_xlim([20, 325])
    ax.set_ylim([0, 30])
    ax.set_ylabel('Mean S(n)')
    ax.set_xlabel('Distance from Centre \n of Tumor (cells)')

    cmap_ax = ax  # for later use

    """
    Large Tumor of Size 10^8 but no new mutations are introduced after 10^6.
    """
    d_append = '_108'

    ax = fig.add_subplot(2, 5, 6 + i)
    plot_it(data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d=d + d_append), ax)
    plot_density(ax, folder, seeds, d=d + d_append)
    ax.set_title('(b) $N=10^8, d=0.' + d[1:] + '$')
    ax.set_xlim([20, 325])
    ax.set_ylim([0, 30])
    ax.set_ylabel('Mean S(n)')


    fig.tight_layout(pad=1, w_pad=0.5)

    subaxes = inset_axes(cmap_ax,
                         width="90%",  # width = 30% of parent_bbox
                         height=0.15,  # height : 1 inch
                         loc=2)

    create_colorbar(subaxes, with_yellow=True)
    return fig


FOLDER_LOCATION = "../model/experiments/u0.01transition"
SEEDS = ['0', '1', '2', '3', '4', '5']
S_plot_transition(FOLDER_LOCATION, SEEDS).savefig('S-plot-transitions.pdf')
