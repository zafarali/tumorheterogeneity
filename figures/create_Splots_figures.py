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

ALL_SEEDS = ['10','100','102','15','3','3318','33181','33185','33186','34201810','342018101','342018102','8','9','99']

sns.set_style('white')
sns.set_context('paper', font_scale=1.5)
from glob import glob
YELLOW = sns.xkcd_rgb["amber"]
params = {'mathtext.default': 'regular' }
plt.rcParams.update(params)
CTC_SIZES = [(2,7), (8,12), (13,17), (18,22),(23,30)]

ALTERNATE_COLORS = sns.xkcd_palette(["vomit green", "orangish brown", "azure", "dark pink"])
small_c, med_c, big_c, biggest_c = ALTERNATE_COLORS



# load data from HCC tumor

HCC = pd.read_csv('./empirical/HCtumordata.csv')
HCC = HCC.sort(columns='r')

# REPLICATE 1
# SELECTED_1_0_0  = '../model/experiments/u0.01875/1_0_0_outs_10/Mar1pipe_out_Fri_Mar__2_20_42_31_2018'
# SELECTED_1_0_02 = '../model/experiments/u0.01875/1_0_02_outs_10/Mar1pipe_out_Fri_Mar__2_21_18_46_2018'

# REPLICATE 2
# SELECTED_1_0_0  = '../model/experiments/u0.01875/1_0_0_outs_99/Mar1pipe_out_Sat_Mar__3_01_01_33_2018'
# SELECTED_1_0_02 = '../model/experiments/u0.01875/1_0_02_outs_99/Mar1pipe_out_Sat_Mar__3_01_44_36_2018'

# VARIABLE
SELECTED_1_0_0 = sys.argv[1]
SELECTED_1_0_02 = sys.argv[2]
SPECIAL_SEED = sys.argv[3]

# S_list_ordered.npy
def empirical_compare_plot(root_folder, seeds):
    NT_folder = root_folder + '/1_0'
    ST_folder = root_folder + '/1_1'
    T_folder = root_folder + '/1_0'

    HCC = pd.read_csv('./empirical/HCtumordata.csv')
    HCC = HCC.sort(columns='r')

    f = plt.figure(figsize=(10, 8))

    ax = f.add_subplot(221)
    ax.set_xlabel('Distance from Centre of Tumor (cm)')
    ax1_r = ax.twinx()
    ax1_r.fill_between(HCC.r, HCC.purity, alpha=0.09, color=(0, 0, 0))
    ax1_r.set_ylabel('Tumor Purity', rotation=-90, labelpad=15)
    x = np.linspace(HCC.r.min(), HCC.r.max(), 43)
    #     y = exp_avg(x)
    dist_sampled_ = sm.add_constant(HCC.r)
    model = sm.OLS(HCC.SNV_corrected, dist_sampled_)
    results = model.fit()

    y = results.predict(sm.add_constant(x))
    ax.plot(x, y, label='Regression (p=' + str(round(results.pvalues[1], 3)) + ')', color=biggest_c)
    ax.scatter(HCC.r, HCC.SNV_corrected, label='Raw Counts', color=biggest_c)
    ax.set_xlim([0, HCC.r.max() + 0.1])
    ax.set_ylim(bottom=0)
    ax.set_ylabel('# Somatic Mutations')
    ax.legend(loc=3)
    ax.set_title('(a) Somatic Mutations in HCC Patient\n')

    folder = root_folder + '/1_0'
    ax = f.add_subplot(223)
    plot_density(ax, folder, seeds, d='0')

    S = np.load(SELECTED_1_0_0+'/S_list01_big.npy')[:,3]
    distances = np.load(SELECTED_1_0_0+'/dist_01_big.npy')[:,3]
    sampled = random.sample(xrange(0, 100), 23)
    s_sampled = S[sampled]
    dist_sampled = distances[sampled]
    ax.scatter(dist_sampled, s_sampled, label='Raw Counts', color=biggest_c)
    #     exp_avg = np.vectorize(exponential_mean_function(dist_sampled, s_sampled ,k=0.001))
    x = np.linspace(distances.min(), distances.max(), 43)
    #     y = exp_avg(x)
    dist_sampled_ = sm.add_constant(dist_sampled)
    model = sm.OLS(s_sampled, dist_sampled_)
    results = model.fit()

    y = results.predict(sm.add_constant(x))
    ax.plot(x, y, label='Regression (p=' + str(round(results.pvalues[1], 3)) + ')', color=biggest_c)
    ax.set_ylabel('# Somatic Mutations')
    ax.set_xlabel('Distance from Centre of Tumor (cells)')
    ax.set_title('(c) No Turnover (>10% frequency)\n')
    ax.set_xlim([0, 325])
    ax.legend(loc=2)
    #     ax.set_ylim([0,60])

    ax = f.add_subplot(224)

    plot_density(ax, folder, seeds, d='02')
    S = np.load(SELECTED_1_0_02+'/S_list01_big.npy')[:,3]
    distances = np.load(SELECTED_1_0_02+'/dist_01_big.npy')[:,3]
    sampled = random.sample(xrange(0, 100), 23)
    s_sampled = S[sampled]
    dist_sampled = distances[sampled]
    ax.scatter(dist_sampled, s_sampled, label='Raw Counts', color=biggest_c)
    #     exp_avg = np.vectorize(exponential_mean_function(dist_sampled, s_sampled ,k=0.001))
    x = np.linspace(distances.min(), distances.max(), 43)
    #     y = exp_avg(x)
    dist_sampled_ = sm.add_constant(dist_sampled)
    model = sm.OLS(s_sampled, dist_sampled_)
    results = model.fit()

    y = results.predict(sm.add_constant(x))
    ax.plot(x, y, label='Regression (p=' + str(round(results.pvalues[1], 3)) + ')', color=biggest_c)
    ax.set_xlabel('Distance from Centre of Tumor (cells)')

    ax.legend(loc=2)

    ax.set_xlim([0, 350])
    ax.set_ylabel('# Somatic Mutations')
    ax.set_title('(d) Turnover (>10% frequency)\n')

    cmapa = create_colormap()

    ax = f.add_subplot(222)
    alternating_power_plot(json.load(open('./turnover_power.json', 'r')), ax)
    ax.set_title('(b) Power Analysis \n (p<0.01)')

    f.tight_layout(h_pad=1.0, w_pad=0.7)
    return f


"""
Figure comparing empirical and actual tumor estimates
"""

empirical_compare_plot('../model/experiments/u0.01875',seeds=SPECIAL_SEED).savefig('fig03.pdf')

"""
Supplementary figure with power analysis for no turnover model
"""
alt_powers = json.load(open('./noturnover_power.json', 'r'))
f = plt.figure()
ax = f.add_subplot(111)
alternating_power_plot(alt_powers, ax)
ax.set_title('No Turnover (f>10%)\n Significance Threshold = 0.01')
ax.set_ylim([0,1.01])
f.savefig('./noturnover_power.pdf')


"""
Figure 2 with multiple S_plots.
"""


def S_plot_paper(root_folder, seeds, k=0.01, pts=100, cutoff='00'):
    fig = plt.figure(figsize=(13, 3))

    ax = fig.add_subplot(141)
    folder = root_folder + '/1_0'
    d2 = '0'

    plot_it(data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d=d2), ax)
    #     plot_it(data_to_plot(folder, seeds, yaxis='S_list'+cutoff+'_big', mode=2, d=d2), ax)
    plot_density(ax, folder, seeds, d=d2)

    ax.set_title('(a) No Turnover')
    ax.set_xlim([50, 325])
    ax.set_ylim([0, 30])
    ax.set_ylabel('Mean S(n)')
    ax.set_xlabel('Distance from Centre \n of Tumor (cells)')

    cmap_ax = ax  # for later use

    d = '005'
    ax = fig.add_subplot(142)
    folder = root_folder + '/1_0'

    plot_it(data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d=d), ax)
    #     plot_it(data_to_plot(folder, seeds, yaxis='S_list'+cutoff+'_big', mode=2, d=d), ax)
    plot_density(ax, folder, seeds, d=d)

    ax.set_title('(b) Turnover, d=0.' + d[1:])
    ax.set_xlim([50, 325])
    ax.set_ylabel('Mean S(n)')
    ax.set_xlabel('Distance from Centre \n of Tumor (cells)')
    #     ax.set_ylim([0,1])

    d = '01'
    ax = fig.add_subplot(143)
    folder = root_folder + '/1_0'
    plot_it(data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d=d), ax)
    #     plot_it(data_to_plot(folder, seeds, yaxis='S_list'+cutoff+'_big', mode=2, d=d), ax)
    plot_density(ax, folder, seeds, d=d)

    ax.set_title('(c) Turnover, d=0.' + d[1:])
    ax.set_ylabel('Mean S(n)')
    x_max = 350 if d == '02' else 325
    x_max = 750 if d == '065' else x_max
    ax.set_xlim([50, x_max])
    ax.set_xlabel('Distance from Centre \n of Tumor (cells)')

    d = '02'
    ax = fig.add_subplot(144)
    folder = root_folder + '/1_0'

    plot_it(data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d=d), ax)
    #     plot_it(data_to_plot(folder, seeds, yaxis='S_list'+cutoff+'_big', mode=2, d=d), ax)
    plot_density(ax, folder, seeds, d=d)

    ax.set_title('(d) Turnover, d=0.' + d[1:])
    ax.set_ylabel('Mean S(n)')
    x_max = 350 if d == '02' else 325
    x_max = 750 if d == '065' else x_max
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


print('---SPLOT PAPER---')
S_plot_paper('../model/experiments/u0.01875/',ALL_SEEDS, cutoff='00').savefig('fig02.pdf')
print('---SPLOT PAPER DONE---')

"""
Supplementary figure for SPlots
"""
def S_plot_paper_supp(root_folder, seeds, k=0.01, pts=100, cutoff='00'):
    fig = plt.figure(figsize=(13, 3))


    d = '005'
    ax = fig.add_subplot(141)

    cmap_ax = ax  # for later use

    folder = root_folder + '/1_1'

    plot_it(data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d=d), ax)
    #     plot_it(data_to_plot(folder, seeds, yaxis='S_list'+cutoff+'_big', mode=2, d=d), ax)
    plot_density(ax, folder, seeds, d=d)

    ax.set_title('(b) Surface Turnover, d=0.' + d[1:])
    ax.set_xlim([50, 325])
    ax.set_ylabel('Mean S(n)')
    ax.set_xlabel('Distance from Centre \n of Tumor (cells)')
    #     ax.set_ylim([0,1])

    d = '01'
    ax = fig.add_subplot(142)
    folder = root_folder + '/1_1'
    plot_it(data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d=d), ax)
    #     plot_it(data_to_plot(folder, seeds, yaxis='S_list'+cutoff+'_big', mode=2, d=d), ax)
    plot_density(ax, folder, seeds, d=d)

    ax.set_title('(c) Surface Turnover, d=0.' + d[1:])
    ax.set_ylabel('Mean S(n)')
    x_max = 350 if d == '02' else 325
    x_max = 750 if d == '065' else x_max
    ax.set_xlim([50, x_max])
    ax.set_xlabel('Distance from Centre \n of Tumor (cells)')

    d = '02'
    ax = fig.add_subplot(143)
    folder = root_folder + '/1_1'

    plot_it(data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d=d), ax)
    #     plot_it(data_to_plot(folder, seeds, yaxis='S_list'+cutoff+'_big', mode=2, d=d), ax)
    plot_density(ax, folder, seeds, d=d)

    ax.set_title('(d) Surface Turnover, d=0.' + d[1:])
    ax.set_ylabel('Mean S(n)')
    x_max = 350 if d == '02' else 325
    x_max = 750 if d == '065' else x_max
    ax.set_xlim([50, x_max])
    ax.set_xlabel('Distance from Centre \n of Tumor (cells)')
    #     ax.set_ylim([0,50])

    d = '065'
    ax = fig.add_subplot(144)
    folder = root_folder + '/1_1'

    plot_it(data_to_plot(folder, seeds, yaxis='S_list_ordered', mode=2, d=d), ax)
    #     plot_it(data_to_plot(folder, seeds, yaxis='S_list'+cutoff+'_big', mode=2, d=d), ax)
    plot_density(ax, folder, seeds, d=d)

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

print '---SPLOT SUPP---'
S_plot_paper_supp('../model/experiments/u0.01875/',ALL_SEEDS, cutoff='00').savefig('Splot_supplementary.pdf')
print('---SPLOT SUPP DONE---')

"""
Supplementary figure showing turnover ad d=0.65
"""
print('----SUPP FIG 3:turnover d=0.65---')
f = plt.figure()
ax = f.add_subplot(111)
plot_it(data_to_plot('../model/experiments/u0.01875/1_0', ALL_SEEDS , yaxis='S_list_ordered', mode=2, d='065'), ax)
plot_density(ax, '../model/experiments/u0.01875/1_0', ALL_SEEDS , d='065')
subaxes = inset_axes(ax,
                    width="90%", # width = 30% of parent_bbox
                    height=0.15, # height : 1 inch
                    loc=2)

create_colorbar(subaxes,with_yellow=True)
ax.set_ylim([0.01, 450])
ax.set_xlim([10, 725])
ax.set_xlabel('Distance from Centre of Tumor (cells)')
ax.set_ylabel('Mean S(n)')
ax.set_title('Turnover d=0.65')
f.savefig('t065_supp.pdf')
print('---SUPP FIG END---')