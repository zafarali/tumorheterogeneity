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


ALTERNATE_COLORS = sns.xkcd_palette(["vomit green", "orangish brown", "azure", "dark pink"])
small_c, med_c, big_c, biggest_c = ALTERNATE_COLORS

"""
Fanning Plots
"""

N = 10**8
b = 0.69

# plots the "fanning_v2" plots
def plot_fanning_v2(root_folder, seeds, k=0.01, pts=100, cutoff='00', d_append='', gamma=0.5, mu=0.02):

    folder = root_folder + '/0_0'


    R_f = lambda d : (3 * N / (4. * np.pi * (1 - d/b)))**(1.0/3.0)


    # turnover value = no-turnover-value + n mu (Rf-R) d/(gamma b)
    fig = plt.figure(figsize=(8, 3))

    ax = fig.add_subplot(111)

    No_turnover = data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d='0'+d_append, max_dist=300)
    Turnover005 = data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d='005'+d_append, max_dist=300)
    Turnover01 = data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d='01'+d_append, max_dist=300)
    Turnover02 = data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d='02'+d_append, max_dist=300)
    Turnover065 = data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d='065'+d_append, max_dist=300)


    x,y = [No_turnover[(1,)][0], No_turnover[(1,)][1] + N*mu*0.05*(R_f(0.05) - Turnover005[(1,)][0])/(gamma * b) ]
    ax.plot(Turnover005[(1,)][0], Turnover005[(1,)][1], '--', color=small_c, label='d=0.05, observed')
    ax.plot(x, y, '-', color=small_c, label='d=0.05, predicted')


    x,y = [No_turnover[(1,)][0], No_turnover[(1,)][1] + N*mu*0.1*(R_f(0.1) - Turnover01[(1,)][0])/(gamma * b) ]
    ax.plot(Turnover01[(1,)][0], Turnover01[(1,)][1], '--', color=med_c, label='d=0.1, observed')
    ax.plot(x, y, '-', color=med_c, label='d=0.1, predicted')



    x,y = [No_turnover[(1,)][0], No_turnover[(1,)][1] + N*mu*0.2*(R_f(0.2) - Turnover02[(1,)][0])/(gamma * b) ]
    ax.plot(Turnover02[(1,)][0], Turnover02[(1,)][1], '--', color=big_c, label='d=0.2, observed')
    ax.plot(x, y, '-', color=big_c, label='d=0.2, predicted')



    x,y = [No_turnover[(1,)][0], No_turnover[(1,)][1] + N*mu*0.65*(R_f(0.65) - Turnover065[(1,)][0])/(gamma * b) ]
    ax.plot(Turnover065[(1,)][0], Turnover065[(1,)][1], '--', color=biggest_c, label='d=0.65, observed')
    ax.plot(x, y, '-', color=biggest_c, label='d=0.65')

    ax.set_xlabel('Distance from COM of Tumor')
    ax.set_ylabel('S(1)')
    ax.legend()
    ax.semilogy(basey=10)
    fig.tight_layout(h_pad=1)
    # return fig
    return fig




plot_fanning_v2('../model/experiments/u0.01/',ALL_SEEDS).savefig('./Splot-radius-increase-fanning0.02.pdf')
plot_fanning_v2('../model/experiments/u0.01/',ALL_SEEDS, mu=0.01).savefig('./Splot-radius-increase-fanning0.01.pdf')




