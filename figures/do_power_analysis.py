import sys
import numpy as np
from sklearn.linear_model import LinearRegression
import random
import statsmodels.api as sm
import json
from glob import glob

def calculate_power_graph_big(distances, S, significance_threshold=0.05):
    """
    Used to do a power analysis for big clusters
    of size 100,1000,10000,20000
    :param distances:
    :param S:
    :param ctc_min_size:
    :param ctc_max_size:
    :param significance_threshold:
    :return:
    """
    total_samples = distances.shape[0]
    significance_count_neg = np.zeros(total_samples)
    significance_count_pos = np.zeros(total_samples)

    for n in range(2, total_samples):
        for i in range(100):
            sample_idx = random.sample(xrange(0, total_samples), n)  # select indexes to look at
            # column 2 contains all the CTCs we want to look at
            #             ctc_size_selection_idx = 2*np.ones(n)
            sample_S = S[sample_idx].reshape(-1, 1)
            sample_distances = distances[sample_idx].reshape(-1, 1)

            # prepare for regression
            sample_distances = sm.add_constant(sample_distances)
            model = sm.OLS(sample_S, sample_distances)
            results = model.fit()
            #         print results.pvalues
            #         print results.summary()
            # do p-value calculations
            p_value = results.pvalues[1]

            if results.params[1] > 0:
                significance_count_pos[n] += int(p_value < significance_threshold)
            else:
                significance_count_neg[n] += int(p_value < significance_threshold)

    significance_count_pos /= 100.0
    significance_count_neg /= 100.0

    return {'pos': significance_count_pos.tolist(), 'neg': significance_count_neg.tolist()}


def calculate_power_graph(distances, S, ctc_min_size=1, ctc_max_size=2, significance_threshold=0.05):
    """
    Used to do power analysis for small clusters between sizes ctc_min_size and ctc_max_size
    :param distances:
    :param S:
    :param ctc_min_size:
    :param ctc_max_size:
    :param significance_threshold:
    :return:
    """
    total_samples = distances.shape[0]
    significance_count_neg = np.zeros(total_samples)
    significance_count_pos = np.zeros(total_samples)

    for n in range(2, total_samples):
        for i in range(100):
            sample_idx = random.sample(xrange(0, total_samples), n)  # select indexes to look at
            print(sample_idx)
            # select the CTC sizes to look at within the range
            # using ctc_min_size-1 because the 0th column corresponds to the ctc of size 1
            ctc_size_selection_idx = np.random.randint(ctc_min_size - 1, ctc_max_size, size=n)
            sample_S = S[sample_idx, ctc_size_selection_idx].reshape(-1, 1)
            sample_distances = distances[sample_idx].reshape(-1, 1)
            print('Sample_S',sample_S.shape)
            print('sample_distances',sample_distances.shape)

            # prepare for regression
            sample_distances = sm.add_constant(sample_distances)
            model = sm.OLS(sample_S, sample_distances)
            results = model.fit()
            # do p-value calculations
            try:
                p_value = results.pvalues[1]

                if results.params[1] > 0:
                    significance_count_pos[n] += int(p_value < significance_threshold)
                else:
                    significance_count_neg[n] += int(p_value < significance_threshold)

            except IndexError as e:
                print('Index error occured...')
                print results.pvalues

    significance_count_pos /= 100.0
    significance_count_neg /= 100.0

    return {'pos': significance_count_pos.tolist(), 'neg': significance_count_neg.tolist()}



CTC_sizes = [(1,1), (2,7), (8,12), (13,17), (18, 22), (23, 30)]
Biopsies = [100, 1000, 10000, 20000]
frequency_thresholds = ['00', '001', '01']

def conduct_power_analysis(root_folder, seed, name_append='', significance_threshold=0.01):
    """
    Does the actual creation of the json to be dumped
    :param root_folder:
    :param seed:
    :param name_append:
    :param significance_threshold:
    :return:
    """

    power_plots = {}
    folder_prepared = root_folder+'_outs_'+seed+'/Mar*/'

    # do power analysis for the small CTCs
    S = np.load(glob(folder_prepared + '/S_list_ordered.npy')[0])
    distances = np.load(glob(folder_prepared + '/deltas_ordered.npy')[0])


    for (ctc_min_size, ctc_max_size) in CTC_sizes:
        power_plots['power_'+str(ctc_min_size)+'_'+str(ctc_max_size)] = calculate_power_graph(distances,S,\
                                                                                     ctc_min_size,ctc_max_size,\
                                                                                     significance_threshold)

    power_plots['x_small'] = range(0,distances.shape[0])

    for frequency_threshold in frequency_thresholds:
        # do power analysis for the big biopsies
        S = np.load(glob(folder_prepared + '/S_list'+frequency_threshold+'_big.npy'))
        distances = np.load(glob(folder_prepared + '/dist_'+frequency_threshold+'_big.npy'))

        print('S.shape', S.shape)
        print('distances.shape', distances.shape)

        for i, biopsy_size in enumerate(Biopsies):
            S_h = S[:, i] # the ith index corresponds to biopsy of size biopsy_size
            distance_h = distances[:, i]
            print('S_h', S_h)
            print('distance_h', distance_h)
            power_plots['power_'+str(biopsy_size)+'_'+str(frequency_threshold)] = calculate_power_graph_big(distance_h, S_h, \
                                                                                                            significance_threshold=significance_threshold)

        power_plots['x_biopsy_'+str(frequency_threshold)] = range(0,distance_h.shape[0])

    return power_plots


root_folder = sys.argv[1] # should be the path to the sims
seed = sys.argv[2] # should be the seed we want to check

# conduct for moderate turnover models:
power_plots = conduct_power_analysis(root_folder+'/1_0_005', seed)
json.dump(power_plots, open('./turnover_power_005.json', 'w'))
power_plots = conduct_power_analysis(root_folder+'/1_0_01', seed)
json.dump(power_plots, open('./turnover_power_01.json', 'w'))
power_plots = conduct_power_analysis(root_folder+'/1_0_02', seed)
json.dump(power_plots, open('./turnover_power_02.json', 'w'))

# conduct for no turnover models
power_plots = conduct_power_analysis(root_folder+'/1_0_0', seed)
json.dump(power_plots, open('./noturnover_power.json', 'w'))
