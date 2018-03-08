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

pi = np.pi
a = 1
R = 10
# this is what we will use to find the roots of to get the radius
# for a specific frequency value
f = np.linspace(0.000001,0.20,num=1000)

r_f = lambda r, f, R, a: f - (a**2/(4*pi))*(1/r**2 - r/R**3)
r_f_ = np.vectorize(r_f)
def pf_3d(f, R, a):
    # find the root
    r_ = optimize.root(r_f, 1, args=(f, R, a)).x[0]
    return ((8*pi**2)/a**2)*(r_**5/(1+r_**3/(2*R**3)))

pf_3d_ = np.vectorize(pf_3d)

fig = plt.figure(figsize=(10, 4))
ax = fig.add_subplot(121)


mappings = [ '../model/experiments/u0.01/0_0_0_outs*']
colors_ = ['gray']
labels_ = ['Neutral Simulation']

fig = plt.figure(figsize=(14,4))
f = np.linspace(0.00001,0.21,num=1000)

ax = fig.add_subplot(131)

freq_plot(ax,
          mappings,
          colors_=colors_,
          labels_=labels_,
          neutral=True)

for R in [1, 5, 10, 15, 20, 50, 100, 200]:
    pdf = pf_3d_(f, R, a)
    not_nans = np.logical_not(np.isnan(np.log10(pdf)))
    ax.plot(np.log10(f)[not_nans], np.log10(pdf)[not_nans], label='R='+str(R))

ax.plot(np.log10(f)[not_nans], np.log10(f**(-2.5)), '--', label='$f^{-2.5}$')
ax.plot(np.log10(f)[not_nans], np.log10(f**(-2)), '--', label='$f^{-2}$')
ax.plot(np.log10(f)[not_nans], np.log10(f**(-1)), '--', label='$f^{-1}$')
ax.set_xlabel('$log_{10}(f)$')
ax.set_ylabel('$log_{10}(p(f))$')
sns.despine()
ax.legend()
ax.set_title('a='+str(a))
R = 20

ax = fig.add_subplot(122)
freq_plot(ax,
          mappings,
          colors_=colors_,
          labels_=labels_,
          neutral=True)
for a in [1, 3, 5, 10]:
    pdf = pf_3d_(f, R, a)
    not_nans = np.logical_not(np.isnan(np.log10(pdf)))
    ax.plot(np.log10(f)[not_nans], np.log10(pdf)[not_nans], label='a='+str(a))
# ax.plot(np.log10(f)[not_nans], np.log10(f**(-2.5)), '--', label='$f^{-2.5}$')
# ax.plot(np.log10(f)[not_nans], np.log10(f**(-2)), '--', label='$f^{-2}$')
# ax.plot(np.log10(f)[not_nans], np.log10(f**(-1)), '--', label='$f^{-1}$')
ax.set_xlabel('$log_{10}(f)$')
ax.set_ylabel('$log_{10}(p(f))$')
sns.despine()
ax.set_title('R='+str(R))
ax.legend()
fig.tight_layout()
fig.savefig('./Analytic.pdf')

