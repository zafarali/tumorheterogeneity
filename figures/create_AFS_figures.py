"""
Creates the Allele Frequency Spectra figures for the paper
Python2.7. Note a lot of code in this
is bad coding practice.
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

"""
Allele Frequency Spectra. 
Figure 1
"""

root_folder = '../model/experiments/u0.01/'
death_rate = '005'
mappings = [ root_folder+'1_0_0_*',
            root_folder+'1_1_'+death_rate+'_*',
           root_folder+'1_0_'+death_rate+'_*']


# f = plt.figure(figsize=(6,9))
f = plt.figure(figsize=(12,5))

ax = f.add_subplot(121)

freq_plot(ax, mappings)
# plt.title('Frequency Spectra',fontsize=12)
ax.set_xlabel('$log_{10}(frequency)$')
ax.set_ylabel('$log_{10}(count)$')
ax.set_xlim([-2.6, 0.05])
ax.legend(fontsize=12,loc=(0.4,0.5))
sns.despine()
ax.set_title('(a) d=0.'+death_rate[1:])

death_rate = '065'

mappings = [ root_folder+'1_0_0_*',
            root_folder+'1_1_'+death_rate+'_*',
           root_folder+'1_0_'+death_rate+'_*']

ax = f.add_subplot(122)
freq_plot(ax, mappings)
ax.set_xlabel('$log_{10}(frequency)$')
ax.set_ylabel('$log_{10}(count)$')
ax.set_xlim([-2.6, 0.05])
ax.legend(fontsize=12,loc=(0.4,0.5))
sns.despine()
ax.set_title('(b) d=0.'+death_rate[1:])
plt.tight_layout(h_pad=1)
plt.savefig('./freqspec001.pdf')


"""
Allele Frequency Spectra. 
Figure S1
"""

death_rate = '01'
mappings = [ root_folder+'1_0_0_*',
            root_folder+'1_1_'+death_rate+'_*',
           root_folder+'1_0_'+death_rate+'_*']


f = plt.figure(figsize=(12,5))
ax = f.add_subplot(121)

freq_plot(ax, mappings)
# plt.title('Frequency Spectra',fontsize=12)
ax.set_xlabel('$log_{10}(frequency)$')
ax.set_ylabel('$log_{10}(count)$')
ax.set_xlim([-2.6, 0.05])
ax.legend(fontsize=12,loc=(0.4,0.5))
sns.despine()
ax.set_title('d=0.'+death_rate[1:]+'\n')
# ax.set_title('d=0')
# ax.savefig('./freqspec'+death_rate+'.pdf')

death_rate = '02'
mappings = [ root_folder+'1_0_0_*',
            root_folder+'1_1_'+death_rate+'_*',
           root_folder+'1_0_'+death_rate+'_*']

ax = f.add_subplot(122)
freq_plot(ax, mappings)
ax.set_xlabel('$log_{10}(frequency)$')
ax.set_ylabel('$log_{10}(count)$')
ax.set_xlim([-2.6, 0.05])
ax.legend(fontsize=12,loc=(0.4,0.5))
sns.despine()
ax.set_title('d=0.'+death_rate[1:]+'\n')
plt.tight_layout(h_pad=1)
# plt.suptitle('Allele')
plt.savefig('./AFS-d001.pdf')