import numpy as np
from sklearn.linear_model import LinearRegression
import random
import statsmodels.api as sm
import json


NOTURNOVER_FOLDER  = '../model/experiments/u0.01875/1_0_0_outs_10/Mar1pipe_out_Fri_Mar__2_20_42_31_2018/'
TURNOVER_FOLDER  = '../model/experiments/u0.01875/1_0_02_outs_10/Mar1pipe_out_Fri_Mar__2_21_18_46_2018/'

FOLDER = TURNOVER_FOLDER

def calculate_power_graph_big(distances, S, ctc_min_size=1, ctc_max_size=2, significance_threshold=0.05):
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


# ctc_min_size, ctc_max_size = 2, 7
# significance_threshold = 0.05

# power_calculation

def calculate_power_graph(distances, S, ctc_min_size=1, ctc_max_size=2, significance_threshold=0.05):
    total_samples = distances.shape[0]
    significance_count_neg = np.zeros(total_samples)
    significance_count_pos = np.zeros(total_samples)

    for n in range(2, total_samples):
        for i in range(100):
            sample_idx = np.random.randint(0, total_samples, size=n)  # select indexes to look at
            # select the CTC sizes to look at within the range
            # using ctc_min_size-1 because the 0th column corresponds to the ctc of size 1
            ctc_size_selection_idx = np.random.randint(ctc_min_size - 1, ctc_max_size, size=n)
            sample_S = S[sample_idx, ctc_size_selection_idx].reshape(-1, 1)
            sample_distances = distances[sample_idx].reshape(-1, 1)

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
                print
                results.pvalues

    significance_count_pos /= 100.0
    significance_count_neg /= 100.0

    return {'pos': significance_count_pos.tolist(), 'neg': significance_count_neg.tolist()}


S = np.load(FOLDER + 'S_list_ordered.npy')
distances = np.load(FOLDER + 'deltas_ordered.npy')[:, 0]



pow2001 = calculate_power_graph(distances,S,1,1,0.01)
pow7001 = calculate_power_graph(distances,S,2,7,0.01)
pow12001 = calculate_power_graph(distances,S,8,12,0.01)
pow17001 = calculate_power_graph(distances,S,13,17,0.01)
pow22001 = calculate_power_graph(distances,S,18,22,0.01)
powbig001 = calculate_power_graph(distances,S,23,30,0.01)
x = range(0,distances.shape[0])

S_h = np.load(FOLDER + 'S_list01_big.npy')[:,3]
distances_h = np.load(FOLDER + '/dist_01_big.npy')[:,3]

powhuge001 = calculate_power_graph_big(distances_h,S_h,significance_threshold=0.01)
x_h = range(0,distances_h.shape[0])

powhuge001_10000 = powhuge001

S_h = np.load(FOLDER + 'S_list01_big.npy')[:,1]
distances_h = np.load(FOLDER + 'dist_01_big.npy')[:,1]

powhuge001_1000 = calculate_power_graph_big(distances_h,S_h,significance_threshold=0.01)

S_h = np.load(FOLDER + 'S_list01_big.npy')[:,0]
distances_h = np.load(FOLDER + 'dist_01_big.npy')[:,0]

powhuge001_100 = calculate_power_graph_big(distances_h,S_h,significance_threshold=0.01)

# temporary use 'pos' for backward compatability
power_plots = {
    'pow2001':pow2001,
    'pow7001':pow7001,
    'pow12001':pow12001,
    'pow17001':pow17001,
    'pow22001':pow22001,
    'powbig001':powbig001,
    'powhuge10000':powhuge001_10000,
    'powhuge1000':powhuge001_1000,
    'powhuge100':powhuge001_100,
    'x':x,
    'x_h':x_h
}

json.dump(power_plots, open('./turnover_power.json', 'w'))



FOLDER = NOTURNOVER_FOLDER

S = np.load(FOLDER + 'S_list_ordered.npy')
distances = np.load(FOLDER + 'deltas_ordered.npy')[:, 0]



pow2001 = calculate_power_graph(distances,S,1,1,0.01)
pow7001 = calculate_power_graph(distances,S,2,7,0.01)
pow12001 = calculate_power_graph(distances,S,8,12,0.01)
pow17001 = calculate_power_graph(distances,S,13,17,0.01)
pow22001 = calculate_power_graph(distances,S,18,22,0.01)
powbig001 = calculate_power_graph(distances,S,23,30,0.01)
x = range(0,distances.shape[0])

S_h = np.load(FOLDER + 'S_list01_big.npy')[:,3]
distances_h = np.load(FOLDER + '/dist_01_big.npy')[:,3]

powhuge001 = calculate_power_graph_big(distances_h,S_h,significance_threshold=0.01)
x_h = range(0,distances_h.shape[0])

powhuge001_10000 = powhuge001

S_h = np.load(FOLDER + 'S_list01_big.npy')[:,1]
distances_h = np.load(FOLDER + 'dist_01_big.npy')[:,1]

powhuge001_1000 = calculate_power_graph_big(distances_h,S_h,significance_threshold=0.01)

S_h = np.load(FOLDER + 'S_list01_big.npy')[:,0]
distances_h = np.load(FOLDER + 'dist_01_big.npy')[:,0]

powhuge001_100 = calculate_power_graph_big(distances_h,S_h,significance_threshold=0.01)


# temporary use 'pos' for backward compatability
power_plots = {
    'pow2001':pow2001,
    'pow7001':pow7001,
    'pow12001':pow12001,
    'pow17001':pow17001,
    'pow22001':pow22001,
    'powbig001':powbig001,
    'powhuge10000':powhuge001_10000,
    'powhuge1000':powhuge001_1000,
    'powhuge100':powhuge001_100,
    'x':x,
    'x_h':x_h
}

json.dump(power_plots, open('./noturnover_power.json', 'w'))