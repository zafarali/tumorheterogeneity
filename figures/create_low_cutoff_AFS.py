import random
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy import optimize
sns.set_style('white')
sns.set_context('paper', font_scale=1.5)
from glob import glob
from figtools import *
from fusco import SFS


pi = np.pi
mu = 0.02
alpha = 30
fusco_alpha = 0.55
fusco_beta = 2.3

mappings = [ '../model/experiments/low_cutoff/0_0_0_low_outs*']
colors_ = ['gray']
labels_ = ['Neutral, $N=10^6, d=0$']

fig = plt.figure()

ax = fig.add_subplot('111')
freq_plot(ax,
          mappings,
          colors_=colors_,
          labels_=labels_,
          neutral=True,
          noS=True,
          calculate_slopes=False,
          cutoff_line_at=0.000001)
fusco_support = np.logspace(-6.0, 0, num=1000)
sfs, x_c_prime = SFS(10**6, mu, fusco_alpha, fusco_beta)
sfs = np.vectorize(sfs)
ax.plot(np.log10(fusco_support), np.log10(sfs(2*fusco_support)), '--r', label='Fusco et al. (With Correction)')
ax.set_xlabel('$log_{10}(f)$')
ax.set_ylabel('$log_{10}(p(f))$')
ax.set_ylim(0)
ax.set_xlim([-6.15, 0.])
sns.despine()

ax.legend(loc='upper right')
# fig.tight_layout()
fig.savefig('./low_cutoff.pdf')

