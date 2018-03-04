"""
Creates all the figures for the paper
Python2.7. Note a lot of code in this
is bad coding practice.

To supply data for this file you need to first run `finescaleAFS.py` as follows:

python finescaleAFS.py turnover close PATH_TO_ALL_SIMS
python finescaleAFS.py turnover notclose PATH_TO_ALL_SIMS
python finescaleAFS.py noturnover close PATH_TO_ALL_SIMS
python finescaleAFS.py notturnover notclose PATH_TO_ALL_SIMS

"""
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


f = plt.figure(figsize=(8,6))

ax = f.add_subplot(221)

S = np.load('./turnover_close_0.npy')
ax.semilogx(np.arange(0,len(S),1)[S>0], S[S>0], 'xb', label='>0%')
S = np.load('./turnover_close_001.npy')
ax.semilogx(np.arange(0,len(S),1)[S>0], S[S>0], 'xr', label='>1%')
S = np.load('./turnover_close_01.npy')
ax.semilogx(np.arange(0,len(S),1)[S>0], S[S>0], 'xg', label='>10%')

ax.set_title('Turnover, center')
ax.set_xlabel('CTC Size')
ax.set_ylabel('S')
ax.set_ylim([0,400])


ax = f.add_subplot(222)
S = np.load('./turnover_notclose_0.npy')
ax.semilogx(np.arange(0,len(S),1)[S>0], S[S>0], 'xb', label='>0%')
S = np.load('./turnover_notclose_001.npy')
ax.semilogx(np.arange(0,len(S),1)[S>0], S[S>0], 'xr', label='>1%')
S = np.load('./turnover_notclose_01.npy')
ax.semilogx(np.arange(0,len(S),1)[S>0], S[S>0], 'xg', label='>10%')
ax.set_title('Turnover, edge')
ax.set_xlabel('CTC Size')
ax.set_ylabel('S')
ax.set_ylim([0,400])
# ax.setlegend()



ax = f.add_subplot(223)
S = np.load('./noturnover_close_0.npy')
ax.semilogx(np.arange(0,len(S),1)[S>0], S[S>0], 'xb', label='>0%')
S = np.load('./noturnover_close_001.npy')
ax.semilogx(np.arange(0,len(S),1)[S>0], S[S>0], 'xr', label='>1%')
S = np.load('./noturnover_close_01.npy')
ax.semilogx(np.arange(0,len(S),1)[S>0], S[S>0], 'xg', label='>10%')

ax.set_title('No Turnover, center')
ax.set_xlabel('CTC Size')
ax.set_ylabel('S')
ax.set_ylim([0,80])
# plt.legend()


ax = f.add_subplot(224)
S = np.load('./noturnover_notclose_0.npy')
ax.semilogx(np.arange(0,len(S),1)[S>0], S[S>0], 'xb', label='>0%')
S = np.load('./noturnover_notclose_001.npy')
ax.semilogx(np.arange(0,len(S),1)[S>0], S[S>0], 'xr', label='>1%')
S = np.load('./noturnover_notclose_01.npy')
ax.semilogx(np.arange(0,len(S),1)[S>0], S[S>0], 'xg', label='>10%')

ax.set_title('No Turnover, edge')
ax.set_xlabel('CTC Size')
ax.set_ylabel('S')
ax.set_ylim([0,80])
ax.legend(loc=2)

f.set_tight_layout({'pad':0.5})

f.savefig('SvsN.pdf')

