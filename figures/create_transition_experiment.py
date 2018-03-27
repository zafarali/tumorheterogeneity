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


def S_plot_transition(root_folder, seeds, d='01', stop_mutants='_106'):
    fig = plt.figure(figsize=(16, 5))
    folder = root_folder + '/0_0_'

    """
    SMall tumor of size 10^6
    """

    tumor_size = '_106' # the size of the tumor

    ax = fig.add_subplot(1, 3, 1)
    plot_it(data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d=d+tumor_size+stop_mutants), ax)
    plot_density(ax, folder, seeds, d=d + tumor_size)
    ax.set_title('(a) $N=10^6, d=0.' + d[1:] + '$')
    ax.set_xlim([20, 325])
    ax.set_ylim([0, 30])
    ax.set_ylabel('Mean S(n)')
    ax.set_xlabel('Distance from Centre \n of Tumor (cells)')

    cmap_ax = ax  # for later use

    """
    Medium Tumor of size 2.4*10^7 
    """


    tumor_size = '_107' # the size of the tumor

    ax = fig.add_subplot(1, 3, 2)
    plot_it(data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d=d+tumor_size+stop_mutants), ax)
    plot_density(ax, folder, seeds, d=d + tumor_size)
    ax.set_title('(a) $N=2.4*10^7, d=0.' + d[1:] + '$')
    ax.set_xlim([20, 325])
    ax.set_ylim([0, 30])
    ax.set_ylabel('Mean S(n)')
    ax.set_xlabel('Distance from Centre \n of Tumor (cells)')

    """
    Large Tumor of Size 10^8 
    """
    tumor_size = '_108'

    ax = fig.add_subplot(1, 3, 3)
    plot_it(data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d=d+tumor_size+stop_mutants), ax)
    plot_density(ax, folder, seeds, d=d+tumor_size+stop_mutants)
    ax.set_title('(b) $N=10^8, d=0.' + d[1:] + '$')
    ax.set_xlim([20, 325])
    ax.set_ylim([0, 30])
    ax.set_ylabel('Mean S(n)')
    ax.set_xlabel('Distance from Centre \n of Tumor (cells)')


    fig.tight_layout(pad=1, w_pad=0.5)

    subaxes = inset_axes(cmap_ax,
                         width="90%",  # width = 30% of parent_bbox
                         height=0.15,  # height : 1 inch
                         loc=2)

    create_colorbar(subaxes, with_yellow=True)
    return fig


FOLDER_LOCATION = "../model/experiments/u0.01transition"

SEEDS = ['6', '7', '8']
"""
With a transition point
"""
S_plot_transition(FOLDER_LOCATION, SEEDS, d='01', stop_mutants='_106').savefig('Splot-Transition-d01-106.pdf')
S_plot_transition(FOLDER_LOCATION, SEEDS, d='0', stop_mutants='_106').savefig('Splot-Transition-d0-106.pdf')
S_plot_transition(FOLDER_LOCATION, SEEDS, d='005', stop_mutants='_106').savefig('Splot-Transition-d005-106.pdf')
S_plot_transition(FOLDER_LOCATION, SEEDS, d='02', stop_mutants='_106').savefig('Splot-Transition-d02-106.pdf')
S_plot_transition(FOLDER_LOCATION, SEEDS, d='065', stop_mutants='_106').savefig('Splot-Transition-d065-106.pdf')

"""
Without a transition point
"""
S_plot_transition(FOLDER_LOCATION, SEEDS, d='01', stop_mutants='_108').savefig('Splot-Transition-d01-108.pdf')
S_plot_transition(FOLDER_LOCATION, SEEDS, d='0', stop_mutants='_108').savefig('Splot-Transition-d0-108.pdf')
S_plot_transition(FOLDER_LOCATION, SEEDS, d='005', stop_mutants='_108').savefig('Splot-Transition-d005-108.pdf')
S_plot_transition(FOLDER_LOCATION, SEEDS, d='02', stop_mutants='_108').savefig('Splot-Transition-d02-108.pdf')
S_plot_transition(FOLDER_LOCATION, SEEDS, d='065', stop_mutants='_108').savefig('Splot-Transition-d065-108.pdf')

