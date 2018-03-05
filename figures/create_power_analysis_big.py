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


def alternating_power_plot(power_plots, ax):
    cmapa = create_colormap()
    HUGE_KEYS_ = ['powhuge100', 'powhuge1000', 'powhuge10000']
    KEYS_ = ['pow_S_000', 'pow_S_001', 'pow_S_010']
    COLORS_ = [cmapa(4), cmapa(9), cmapa(14)]
    LABELS_ = ['f>0%', 'f>1%', 'f>10%']

    for k, color, label in zip(KEYS_, COLORS_, LABELS_):

        to_plot = np.array(power_plots[k]['scount'])

        if np.any(to_plot):
            ax.scatter(np.arange(0, to_plot.shape[0], 1), to_plot, color=color, marker='o', label=label)


    ax.legend(loc=1, labelspacing=0, handlelength=0, frameon=True)

    ax.set_xlim([0, 1])
    ax.set_xlabel('Number of samples')
    ax.set_ylabel('Proportion of \nSignificant Regressions')

    return ax

power_plots = json.load(open('./big_powerplot.json', 'r'))

f = plt.figure()
ax = f.add_subplot(111)
ax = alternating_power_plot(power_plots, ax)
ax.set_title('Turnover d=0.02\nSignificance threshold = 0.01')
f.savefig('power_analysis_big.pdf')

