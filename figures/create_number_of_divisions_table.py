"""
Table S1 data comes from this file.
"""
import numpy as np
from glob import glob

mutation_rate = 0.0375  # this is the mutation rate in the simulations

root_folder = '../model/experiments/u0.01875/'

for death_rate in ['005', '01', '02', '065']:

    mappings = [root_folder + '1_0_0_*',
                root_folder + '1_1_' + death_rate + '_*',
                root_folder + '1_0_' + death_rate + '_*']

    print 'DEATH RATE:' + death_rate
    for mapping, model_name in zip(mappings, ['No Turnover', 'Surface Turnover', 'Turnover']):

        replicates = glob(mapping)

        all_data = []
        for replicate_folder in replicates:
            try:
                datafile = glob(replicate_folder + '/*/statistics.csv')[0]
                mean_SNPs_in_clusters = np.mean(pd.read_csv(datafile)['mean_SNPs']) / mutation_rate
                all_data.append(mean_SNPs_in_clusters)
            except Exception as e:
                pass
        #             print 'Exception: '+str(e)
        #             print 'Folder Skipped: '+replicate_folder

        print model_name + ' mean number of divisions: ' + str(np.mean(all_data))
        print model_name + ' SD number of divisions: ' + str(np.std(all_data))

