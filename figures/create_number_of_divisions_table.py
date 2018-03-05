"""
Table S1 data comes from this file.
"""
import numpy as np
from glob import glob
import pandas as pd

mutation_rate = 0.02  # this is the mutation rate in the simulations

root_folder = '../model/experiments/u0.01/'

print '-----TABLE START -----'
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
                datafile = glob(replicate_folder+'/*/S_list_ordered.npy')[0]
                SNPs_per_cell = np.load(datafile).mean(axis=0)[0]
                number_of_divisions = SNPs_per_cell/mutation_rate
                all_data.append(number_of_divisions)
            except Exception as e:
                print('Exception occured:'+str(e))
                pass

        print model_name + ' mean number of divisions: ' + str(np.mean(all_data))
        print model_name + ' SD number of divisions: ' + str(np.std(all_data))

print '----TABLE END ----'