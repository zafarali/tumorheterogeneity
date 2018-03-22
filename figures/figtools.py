"""
Most of the code in this file is just a copy-paste from
different ipython notebooks. It runs in python2.7
"""
import random
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from mpl_toolkits.axes_grid.inset_locator import inset_axes
sns.set_style('white')
sns.set_context('paper', font_scale=1.5)
from glob import glob
import random
import json
import statsmodels.api as sm
YELLOW = sns.xkcd_rgb["amber"]
params = {'mathtext.default': 'regular' }
plt.rcParams.update(params)
CTC_SIZES = [(2,7), (8,12), (13,17), (18,22),(23,30)]

ALTERNATE_COLORS = sns.xkcd_palette(["vomit green", "orangish brown", "azure", "dark pink"])
small_c, med_c, big_c, biggest_c = ALTERNATE_COLORS


# creates a color map that you can use like this: create_colormap(y) --> rgba color value
def create_colormap(colormap_name='cubehelix_r', with_yellow=False, n_colors=8, use_range=(1, 3, 4, 6, 7, 7),
                    start_bin=2, each_bin=5, return_colormap=False):
    colors = sns.color_palette(colormap_name, n_colors=n_colors)
    bin_ = start_bin
    color_mapping = {}
    for color in use_range:  # map through all colors
        for i in xrange(each_bin):
            color_mapping[bin_ + i] = colors[color]
        bin_ = bin_ + each_bin

    if return_colormap:
        col_list = [colors[i] for i in use_range]
        if with_yellow:
            col_list = [YELLOW] + col_list
        #             col_list = [ YELLOW ] + col_list + ALTERNATE_COLORS
        return mpl.colors.ListedColormap(col_list)

    def c_mapper(i):
        print i
        try:
            if i == 1 and with_yellow:
                return YELLOW
            elif i == 100:
                return small_c
            elif i == 1000:
                return med_c
            elif i == 10000:
                return big_c
            elif i == 20000:
                return biggest_c
            else:
                return color_mapping[i]
        except Exception as e:
            print('error')
            return color_mapping[sorted(color_mapping.keys())[-1]]

    return c_mapper

def create_colorbar(subaxes, with_yellow=False):
    cb = mpl.colorbar.ColorbarBase(subaxes, cmap=create_colormap(return_colormap=True, with_yellow=with_yellow),
    #                                 tick=['2-7', '8-13', '14-19', '20-25', '26-30'],
    #                                 extend='both',
                                    # Make the length of each extension
                                    # the same as the length of the
                                    # interior colors:
    #                                 extendfrac='auto',
    #                                 ticks=bounds,
                                    spacing='uniform',
                                    orientation='horizontal')
    cb.set_label('Cluster Sizes')
    cb.ax.get_xaxis().labelpad = 0

    ticklocs = [0.03, 0.2, 0.33, 0.5, 0.65,0.85] # location of ticks
    ticklabels = ['1', '2-7', '8-12', '13-17', '18-22', '23-30', '10^2', '10^3', '10^4', '2*10^4'] # tick labels

    #ticklocs = [0.03, 0.1, 0.21, 0.33, 0.4, 0.5, 0.65, 0.75, 0.85, 0.9] # location of ticks
    #ticklabels = ['1', '2-7', '8-12', '13-17', '18-22', '23-30', '10^2', '10^3', '10^4', '2*10^4'] # tick labels

    if not with_yellow:
        ticklabels = ticklabels[1:]
        ticklocs = list(np.linspace(0.05,0.6,4))+[0.8]

    cb.set_ticks(ticklocs)
    cb.ax.set_xticklabels(ticklabels, fontsize = 9, rotation=30)


def exponential_mean_function(x_i, y_i, k=1):
    """
        @params:
            x_i : all x_i points available
            y_i : the y_i values to average over
        @returns
            a function that can be evaluated for the mean at any point x

    """

    def exponential_mean_at(x):
        """
            @params:
                x: a np.array or float where the average value needs to be returned
            @returns
                <y(x)> : average value at position x
        """
        w = np.exp(-k * (x - x_i) ** 2)
        return np.dot(w, y_i) / np.sum(w)

    return exponential_mean_at


def prepare_data2(folder, yaxis='S_list_ordered', k=0.01, pts=100, max_dist=None):
    """
        Load data from .npy files
    """

    assert yaxis in ['S_list01_big', 'S_list001_big', 'S_list00_big', 'S_list_simulated', 'S_list_ordered',
                     'D_list_ordered', 'cdrv_list_ordered',
                     'csnp_list_ordered'], 'yaxis is not recognized for prepare_data2'

    ISBIG = yaxis in ['S_list01_big', 'S_list001_big', 'S_list00_big']
    dist_file = 'dist_' + yaxis.split('_list')[1] if ISBIG else 'deltas_ordered'
    if len(glob(folder + 'Mar*/' + dist_file + '.npy')) == 0:
        raise Exception('No files found in ' + folder)

    distances = np.load(glob(folder + 'Mar*/' + dist_file + '.npy')[0])

    if not ISBIG:
        distances = distances[:, 0]

    S = np.load(glob(folder + 'Mar*/' + yaxis + '.npy')[0])
    prepped = {}
    all_sizes = [1, 2, 3, 4] if yaxis in ['S_list01_big', 'S_list001_big', 'S_list00_big'] else list([1] + CTC_SIZES)
    for ctc_size in all_sizes:
        slice_range = None
        x = None
        y = None
        color = None
        if type(ctc_size) is int:  # check if this is an int so that we can handle different compared to a range...
            x = distances
            y = S[:, ctc_size - 1]
            if ISBIG:
                y = y[:100]
                x = distances[:100, ctc_size - 1]
            ctc_size = (ctc_size,)
        else:
            slice_range = range(ctc_size[0] - 1, ctc_size[1])
            x = np.repeat(distances, ctc_size[1] - ctc_size[0] + 1)
            y = S[:, slice_range].reshape(1, -1)[0]

        minp = 0
        maxp = distances.max() if max_dist is None else max_dist # define the maximum tumor size
        x_av = np.linspace(minp, maxp, pts)
        moving_avg_function = np.vectorize(exponential_mean_function(x, y, k=k))  # returns a function
        y_av = moving_avg_function(x_av)

        if ISBIG:
            #             print 'idx=',ctc_size[0]
            _ctc_size = 10 ** (ctc_size[0] + 1) if ctc_size[0] != 4 else 2 * 10 ** 4
            #             print 'powered=',_ctc_size
            ctc_size = (_ctc_size,)
        prepped[ctc_size] = [x_av, y_av]

    return prepped


def cross_sim_average(samples):
    to_be_returned = {}
    ctc_sizes = samples[0].keys()
    for ctc_size in ctc_sizes:
        to_be_avged_x = []
        to_be_avged_y = []
        # go through all samples and extract that
        for sample in samples:
            to_be_avged_x.append(sample[ctc_size][0])
            to_be_avged_y.append(sample[ctc_size][1])

        to_be_avged_x = np.array(to_be_avged_x)
        to_be_avged_y = np.array(to_be_avged_y)

        x_avgd = np.mean(to_be_avged_x, axis=0)
        y_avgd = np.mean(to_be_avged_y, axis=0)
        to_be_returned[ctc_size] = [x_avgd, y_avgd]

    return to_be_returned


def data_to_plot(folder, seeds, yaxis='unique_combos', k=0.01, pts=100, mode=2, d='005', max_dist=None):
    """
    :param folder: The folder to look at
    :param seeds:  The seeds to consider
    :param yaxis: The y-axis to plot
    :param k: the smoothing for the kernel
    :param pts: The number of points to return
    :param mode: method to load data (should be 2)
    :param d: the death rate
    :param max_dist: the maximum distance to return data for
    :return:
    """
    samples = [] # all samples will be held here
    # first prep the data for each seed
#     print 'DATA_TO_PLOT - FOLDER:'+folder+' SEEDS:'+str(seeds)+' d:'+d

    for seed in seeds:
        this_folder = folder+'_'+d+'_outs_'+seed+'/'
#         print this_folder
        try:
            if mode == 1:
                raise DeprecationWarning('Cannot use prepare_data1 anymore.')
            else:
                samples.append(prepare_data2(this_folder, yaxis=yaxis, k=k, pts=pts, max_dist=max_dist))
        except Exception as e:
            # TODO: there is an error going on here
            print 'data_to_plot Exception:'+str(e)
            print 'Skipped Folder:'+this_folder
#             raise e
    # produce the cross sample average
    cross_sample_averaged = cross_sim_average(samples)
    return cross_sample_averaged


def plot_it(cross_sample_averaged, ax, linestyle='-'):
    colors = create_colormap(with_yellow=True)
    for ctc_size in cross_sample_averaged.keys():
        x, y = cross_sample_averaged[ctc_size]
        ax.plot(x, y, color=colors(ctc_size[0]), linestyle=linestyle)

        if ctc_size[0] == 1:
            print 'mean(s(1))', np.mean(y)
            print 'std(s(1))', np.std(y)
            print 's(1)_core', np.mean(y[:10])
            print 's(1)_edge', np.mean(y[10:])



def plot_density(ax, folder, seeds, d='005'):
    r_samples = [] # all samples will be held here
    rho_samples = [] # all samples will be held here
    # first prep the data for each seed
#     print 'PLOT DENSITY folder:'+folder+' seeds:'+str(seeds)+' d='+d
    for seed in seeds:
        this_folder = folder+'_'+d+'_outs_'+seed+'/'
#         print this_folder
        try:
            r_samples.append(np.load(glob(this_folder+'/Mar*/r_meaned.npy')[0]))
            rho_samples.append(np.load(glob(this_folder+'/Mar*/rho2.npy')[0]))
        except Exception as e:
            print 'Exception: '+str(e)
            print 'Skipped Folder: '+this_folder
    print len(r_samples)
    print len(rho_samples)
    r = np.mean(np.array(r_samples), axis=0)
    rho = np.mean(np.array(rho_samples), axis=0)

    ax1_r = ax.twinx()
    ax1_r.fill_between(r, rho,alpha=0.09,color=(0,0,0))
    ax1_r.set_ylim([0,1.1])
    ax1_r.set_ylabel('Tumor Cell Density', rotation=-90,labelpad=15)


RED =  sns.xkcd_rgb["pale red"]
# GREEN = sns.xkcd_rgb["medium green"]
# BLUE = sns.xkcd_rgb["denim blue"]

def widths(arr):
    index = 0
    while index < len(arr)-1:
        yield arr[index+1] - arr[index]
        index += 1

def neighbour_iterator(arr):
    index = 0
    while index < len(arr)-1:
        yield (arr[index], arr[index+1])
        index += 1

numbins = 50

binrange = np.logspace(np.log10(0.0001), np.log10(1), num=numbins)

# _, GREEN, RED, _, _, BLUE = sns.color_palette('colorblind')
RED, BLUE, GREEN = sns.xkcd_palette(["amber", "dusty purple", "faded green"])
# RED, BLUE, GREEN = sns.color_palette("cubehelix", 3)
sns.set_context('paper', font_scale=1.5)

# TODO: neutral simultations do not have drv files, need to fork this into two different things
def freq_plot(ax, mappings,
              colors_ = [RED, GREEN, BLUE],
              labels_ = ['No Turnover', 'Surface Turnover', 'Turnover'],
              neutral=False,
              calculate_slopes=False,
              noS=False,
              slope_start=-2.5,
              slope_end=0):
    lines = []
    labels = []
    for mapping, color, model_name in zip(mappings, colors_,labels_):

        replicates = glob(mapping)

        all_xs = []
        all_muts = []
        all_drvs = []

        all_muts_or = []
        all_drvs_or = []

        for replicate_folder in replicates:
            try:
                datafile = glob(replicate_folder + '/all_PMs*.dat')[0]
                muts = pd.read_csv(datafile, sep=' ', names=['Num', 'SNP', 'abundancy'])
                if not neutral: datafile = glob(replicate_folder + '/drv_PMs_*.dat')[0]

                if not neutral:
                    drvs = pd.read_csv(datafile, sep=' ', names=['Num', 'SNP', 'abundancy'])

                #                 print(np.min(muts.abundancy))
                y, x = np.histogram(muts.abundancy, bins=binrange)

                all_muts_or.append(y.copy())
                y = [y_ / float(w_) for y_, w_ in zip(y.tolist(), list(widths(x)))]

                x_meaned = np.mean(np.array(list(neighbour_iterator(x))), axis=1)
                all_muts.append(y)
                all_xs.append(x_meaned)

                if not neutral:
                    y, x = np.histogram(drvs.abundancy, bins=binrange)
                    all_drvs_or.append(y.copy())
                    y = [y_ / float(w_) for y_, w_ in zip(y.tolist(), list(widths(x)))]
                    all_drvs.append(y)

            except Exception as e:
                print 'Could not find file:' + str(replicate_folder) + ' ' + str(e)

        y = np.mean(np.array(all_muts), axis=0)
        y2 = np.mean(np.array(all_drvs), axis=0)
        y_or = np.mean(np.array(all_muts_or), axis=0)
        if not neutral: y2_or = np.mean(np.array(all_drvs_or), axis=0)
        x_meaned = np.mean(np.array(all_xs), axis=0)

        y_keep = y != 0
        if not neutral: y2_keep = y2 != 0

        y1_x_plt = np.log10(x_meaned[y_keep])
        y1_plt = np.log10(y[y_keep]).astype(np.float)

        if not neutral: y2_x_plt = np.log10(x_meaned[y2_keep])
        if not neutral: y2_plt = np.log10(y2[y2_keep]).astype(np.float)

        print model_name + ':'

        if calculate_slopes:
            # calcualte regressions, can probably save some computation by
            # removing log10s everywhere, but ive kept them here for explicitness
            truncation_point = np.argmax(np.isinf(np.log10(y)))
            print(truncation_point)
            x_meaned_truncated_base_10 = np.log10(x_meaned[:truncation_point])
            y_truncated_base_10 = np.log10(y[:truncation_point])
            print(x_meaned_truncated_base_10)
            print(y_truncated_base_10)
            # prepare for regression
            lt25 = x_meaned_truncated_base_10 < slope_start
            x_meaned_lt25 = sm.add_constant(x_meaned_truncated_base_10[lt25])

            model = sm.OLS(y_truncated_base_10[lt25], x_meaned_lt25)
            results = model.fit()
            print('x-values:  10^-4 to 10^'+str(slope_start))
            print('regression of passengers, coefficients', results.params)
            print('regression of passengers, p-values', results.pvalues)
            print('allvalues:', y_truncated_base_10[lt25], x_meaned_lt25)

            gt25 = x_meaned_truncated_base_10 >= slope_start
            x_meaned_selected = x_meaned_truncated_base_10[gt25]
            ltse = x_meaned_selected < slope_end
            x_meaned_selected = x_meaned_selected[ltse]
            x_meaned_gt25 = sm.add_constant(x_meaned_selected)

            model = sm.OLS(y_truncated_base_10[gt25][ltse], x_meaned_gt25)
            results = model.fit()
            print('x-values: 10^'+str(slope_start)+' to 10^'+str(x_meaned_selected[-1]))
            print('regression of passengers, coefficients', results.params)
            print('regression of passengers, p-values', results.pvalues)
            print('allvalues:', y_truncated_base_10[gt25][ltse], x_meaned_gt25)

            if not neutral:
                y2_x_plt_ = sm.add_constant(y2_x_plt)
                model = sm.OLS(y2_plt, y2_x_plt_)
                results = model.fit()
                print('regression on drivers, coefficients', results.params)
                print('regression on drivers, p-values',results.pvalues)

        f_sum = np.around(np.sum(y_or), decimals=2)
        if not neutral: d_sum = np.around(np.sum(y2_or), decimals=2)
        print 'FREQ_SUM:' + str(f_sum)
        if not neutral: print 'DRV_SUM:' + str(d_sum)
        if not neutral: print 'PROP DRV:' + str(d_sum / f_sum)

        label_ = str(model_name)
        if not noS: label_+= ' $S$=' + str(f_sum)
        lines += ax.plot(y1_x_plt, y1_plt, 'o', color=color, alpha=1, label=label_)
        labels += [str(model_name) + ' $S$=' + str(f_sum)]
        if not neutral:
            lines += ax.plot(y2_x_plt, y2_plt, '^', color=color, alpha=0.8,
                             label=str(model_name) + ' (Drivers) $S_d$=' + str(d_sum), markersize=8)
            labels += [str(model_name) + ' (Drivers) $S_d$=' + str(d_sum)]
        ax.plot(np.log10(np.ones(40) * 0.0001), np.linspace(0, 10, num=40), ':', color='gray', linewidth=0.5)


    return lines, labels



CTC_sizes = [(1,1), (2,7), (8,12), (13,17), (18, 22), (23, 30)]


def alternating_power_plot(power_plots, ax, frequency_threshold='00', legend=True):
    cmapa = create_colormap()
    HUGE_KEYS_ = [ 'power_'+s+'_'+frequency_threshold  for s in ['100', '1000', '10000', '20000']]
    CTC_KEYS_ = [ 'power_'+str(ctc_min_size)+'_'+str(ctc_max_size) for (ctc_min_size, ctc_max_size) in CTC_sizes ]
    ALL_KEYS = CTC_KEYS_ + HUGE_KEYS_
    COLORS_ = [YELLOW, cmapa(4), cmapa(9), cmapa(14), cmapa(20), cmapa(24), small_c, med_c, big_c, biggest_c]
    LABELS_ = ['1', '(2,7)', '(8,12)', '(13,17)', '(18,22)', '(23,30)', '100',
               '1000', '10000', '20000']

    for k, color, label in zip(ALL_KEYS, COLORS_, LABELS_):
        x_ = power_plots['x_biopsy_'+frequency_threshold] if k in HUGE_KEYS_ else power_plots['x_small']

        pos = np.array(power_plots[k]['pos'])
        neg = np.array(power_plots[k]['neg'])
        pos_tot = np.sum(pos)
        neg_tot = np.sum(neg)
        marker = 'o' if k in HUGE_KEYS_ else '^'
        if np.any(power_plots[k]['pos']) and pos_tot > neg_tot:
            ax.scatter(x_, power_plots[k]['pos'], color=color, marker=marker, label=label)
        if np.any(power_plots[k]['neg']) and neg_tot > pos_tot:
            ax.scatter(x_, -np.array(power_plots[k]['neg']), color=color, marker=marker, label=label)

    if legend: ax.legend(loc=1, labelspacing=0, handlelength=0, frameon=True)
    ax.plot(np.ones(40) * 23, np.linspace(-2, 2, 40), '--', color='gray')
    ax.set_xlim([0, 60])
    ax.set_ylim([-1.1, 1.1])
    ax.set_xlabel('Number of samples')
    ax.set_ylabel('Signed Proportion of \nSignificant Regressions')

    return ax