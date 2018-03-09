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


RED, BLUE, GREEN = sns.xkcd_palette(["amber", "dusty purple", "faded green"])
sns.set_context('paper', font_scale=1.5)

pi = np.pi
mu = 0.02
alpha = 40

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
ax.set_xlim([-4.1, 0.05])
sns.despine()
ax.set_title('(a) d=0.'+death_rate[1:])
ax.set_ylim(0)

# analytic line
freq_support = np.linspace(0.000001,0.21,num=1000)
ax.plot(np.log10(freq_support), np.log10((alpha*mu/(4*np.sqrt(np.pi)))*freq_support**(-2.5)), '--', label='Analytic Result $alpha=$ '+str(alpha))
ax.legend(fontsize=12,loc=(0.3,0.5))
death_rate = '065'

mappings = [ root_folder+'1_0_0_*',
            root_folder+'1_1_'+death_rate+'_*',
           root_folder+'1_0_'+death_rate+'_*']

ax = f.add_subplot(122)
freq_plot(ax, mappings)
ax.set_xlabel('$log_{10}(frequency)$')
ax.set_ylabel('$log_{10}(count)$')
ax.set_xlim([-4.1, 0.05])
ax.set_ylim(0)
ax.legend(fontsize=12,loc=(0.3,0.6))
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
ax.set_xlim([-4.1, 0.05])
ax.legend(fontsize=10,loc=(0.4,0.5))
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
ax.set_xlim([-4.1, 0.05])
ax.legend(fontsize=10,loc=(0.4,0.5))
sns.despine()
ax.set_title('d=0.'+death_rate[1:]+'\n')
plt.tight_layout(h_pad=1)
# plt.suptitle('Allele')
plt.savefig('./AFS-d001.pdf')


"""
how does death rate in no selection simulations change the from the theoretical prediction?
"""

cmap = create_colormap()
COLORS = [cmap(1), cmap(4), cmap(9), cmap(14), cmap(19)]
sns.set_context('paper', font_scale=1.5)

pi = np.pi
mu = 0.02

root_folder = '../model/experiments/u0.01/'
mappings = [ root_folder+'/0_0_0_outs*',
             root_folder+'/0_0_005_outs*',
             root_folder+'/0_0_01_outs*',
             root_folder+'/0_0_02_outs*',
             root_folder+'/0_0_065_outs*']

fig = plt.figure(figsize=(14,4))
f = np.linspace(0.000001,0.21,num=1000)

ax = fig.add_subplot(131)

freq_plot(ax,
          mappings,
          colors_=COLORS,labels_=['d=0', 'd=0.05', 'd=0.1', 'd=0.2', 'd=0.65'],
          neutral=True)

ax.set_xlabel('$log_{10}(frequency)$')
ax.set_ylabel('$log_{10}($count density)')

ax.plot(np.log10(freq_support), np.log10((mu/(4*np.sqrt(np.pi)))*freq_support**(-2.5)), '--', label='Assuming $E[i^{3/2}]=1$')
ax.plot(np.log10(freq_support), np.log10((15*mu/(4*np.sqrt(np.pi)))*freq_support**(-2.5)), '--', label='Assuming $E[i^{3/2}]=15$', color=BLUE)
ax.plot(np.log10(freq_support), np.log10(freq_support**(-2)), '--', label='$f^{-2}$', color='r')


ax.set_xlim([-4.15, 0.])
ax.set_ylim(bottom=0)
ax.legend(fontsize=8,loc=(0.45,0.3))
sns.despine()
ax.set_title('(a) No Selection')

mappings = [ root_folder+'1_0_0_*',
            root_folder+'1_0_005_*',
           root_folder+'1_0_01_*',
           root_folder+'1_0_02_*',
           root_folder+'1_0_065_*']

ax = fig.add_subplot(132)

freq_plot(ax, mappings, colors_=COLORS,
          labels_=['d=0', 'd=0.05', 'd=0.1', 'd=0.2', 'd=0.65'],neutral=True)

ax.set_xlabel('$log_{10}(frequency)$')
ax.set_ylabel('$log_{10}($count density)')

ax.plot(np.log10(freq_support), np.log10((mu/(4*np.sqrt(np.pi)))*freq_support**(-2.5)), '--', label='Assuming $E[i^{3/2}]=1$')
ax.plot(np.log10(freq_support), np.log10((15*mu/(4*np.sqrt(np.pi)))*freq_support**(-2.5)), '--', label='Assuming $E[i^{3/2}]=15$', color=BLUE)
ax.plot(np.log10(freq_support), np.log10(freq_support**(-2)), '--', label='$f^-2$', color='r')

# ax.plot(np.log10(f), np.log10(f**(-2)*1/(4*np.sqrt(np.pi))), '--', label='f^-2$')
# ax.plot(np.log10(f), np.log10(1/f*(f+1)), '--', label='1/f(f+1)$')

ax.set_xlim([-4.15, 0.])
ax.set_ylim(bottom=0)
ax.legend(fontsize=8, loc=(0.45,0.3))
sns.despine()
ax.set_title('(b) Selection = 1%')


mappings = [ root_folder+'/10_0_0_outs*',
             root_folder+'/10_0_005_outs*',
             root_folder+'/10_0_01_outs*',
             root_folder+'/10_0_02_outs*',
             root_folder+'/10_0_065_outs*']

ax = fig.add_subplot(133)

freq_plot(ax, mappings, neutral=True,
          colors_=COLORS,
          labels_=['d=0', 'd=0.05', 'd=0.1', 'd=0.2', 'd=0.65'])

ax.set_xlabel('$log_{10}(frequency)$')
ax.set_ylabel('$log_{10}($count density)')
freq_support = np.linspace(0.000001,0.21,num=1000)

ax.plot(np.log10(freq_support), np.log10((mu/(4*np.sqrt(np.pi)))*freq_support**(-2.5)), '--', label='Assuming $E[i^{3/2}]=1$')
ax.plot(np.log10(freq_support), np.log10((15*mu/(4*np.sqrt(np.pi)))*freq_support**(-2.5)), '--', label='Assuming $E[i^{3/2}]=15$', color=BLUE)
ax.plot(np.log10(freq_support), np.log10(freq_support**(-2)), '--', label='$f^-2$', color='r')

ax.set_xlim([-4.15, 0.])
ax.set_ylim(bottom=0)
ax.legend(fontsize=8, loc=(0.45,0.3))
sns.despine()
ax.set_title('(c) Selection = 10%')
fig.tight_layout()

# ax.set_title('No Turnover')
plt.savefig('./afs-different_selection_rates.pdf', dpi=300)
