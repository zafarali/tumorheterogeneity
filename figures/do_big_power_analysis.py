import numpy as np
from sklearn.linear_model import LinearRegression
import random
import statsmodels.api as sm
import json


NOTURNOVER_FOLDER  = '../model/experiments/u0.01875/1_0_0_outs_10/Mar1pipe_out_Fri_Mar__2_20_42_31_2018/'
TURNOVER_FOLDER  = '../model/experiments/u0.01875/1_0_02_outs_10/Mar1pipe_out_Fri_Mar__2_21_18_46_2018/'

FOLDER = TURNOVER_FOLDER

def calculate_power_graph_big(distances, S, significance_threshold=0.05):
    total_samples = distances.shape[0]
    significance_count = np.zeros(total_samples)

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

            significance_count[n] += int(p_value < significance_threshold)

    significance_count /= 100.0

    return {'scount': significance_count.tolist()}



S_000 = np.load(FOLDER + 'S_list00_big.npy')[:,3]
S_001 = np.load(FOLDER + 'S_list001_big.npy')[:,3]
S_010 = np.load(FOLDER + 'S_list01_big.npy')[:,3]

distances_h_000 = np.load(FOLDER + '/dist_00_big.npy')[:,3]
distances_h_001 = np.load(FOLDER + '/dist_001_big.npy')[:,3]
distances_h_010 = np.load(FOLDER + '/dist_01_big.npy')[:,3]

pow_S_000 = calculate_power_graph_big(distances_h_000,S_000,significance_threshold=0.01)
pow_S_001 = calculate_power_graph_big(distances_h_001,S_001,significance_threshold=0.01)
pow_S_010 = calculate_power_graph_big(distances_h_010,S_010,significance_threshold=0.01)

power_plots = {
    'pow_S_000':pow_S_000,
    'pow_S_001':pow_S_001,
    'pow_S_010':pow_S_010,
}
json.dump(power_plots, open('./big_powerplot.json', 'w'))
