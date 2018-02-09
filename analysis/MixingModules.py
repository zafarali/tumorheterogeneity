# modules for pipeline
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import csv
import os
import itertools
from BaseObjects import Tumor
from Samplers import SphericalSampler, KDTSphericalSampler
import Statistics
import random
import time
import json
import multiprocessing.dummy as multiprocessing

CPU_COUNT = 12
EPS = np.finfo(float).eps
MAX_CELLS_RANGE = 50
INTERVAL = 5
LARGE_CLUSTERS = [100, 1000, 10000, 20000]
SMALL_CLUSTERS = [2, 6, 27, 50]
CLUSTER_SIZES = SMALL_CLUSTERS + LARGE_CLUSTERS

def load_snps(datafile):
    mutations = pd.read_csv(datafile, sep=' ', names=['Num', 'SNP', 'abundancy'])
    mutations['abundancy'] = mutations['abundancy'] + EPS 
    mutations['sampling_prob'] = mutations['abundancy']/mutations['abundancy'].sum()
    return mutations 

def prepare_tumor_for_mixing_analysis(pipeline):
    """
    Prepares a tumor for the mixing analysis
    """
    pipeline.tumor._make_snp_lookup_efficient()

def perform_mixing_analysis(pipeline):
    """
    Performs mixing analysis. 
    Goes through all snps and
    randomly samples a a set of snps
    so that we can perform mixing analysis.
    """
    pipeline.print2('Performing mixing analysis.')
    # 1) load snps
    mutations = load_snps(pipeline.FILES['all_PMs'])

    pipeline.print2('Loaded all {} mutations. Now sampling'.format(mutations.shape[0]))
    # 2) sample snps according to frequency
    mutation_idx = np.random.choice(mutations['Num'], size=50, replace=False, p=mutations['sampling_prob'])
    mutations_to_analyze = mutations.iloc[mutation_idx]

    pipeline.print2('Sampled 50 mutations. Investigating them individually.')
    print(mutations_to_analyze)
    # print('Sampled SNPS: {}'.format(mutations_to_analyze))

    # 3) perform analysis for each of those snps

    all_results = []
    pool = multiprocessing.Pool(CPU_COUNT)

    count = 0
    pipeline.print2('Starting loop.')
    for row in mutations_to_analyze.iterrows():
        try:
            data = row[1]
            abundancy, snp_id = data.abundancy, int(data.SNP)
            results = _snp_mixing_analysis(pipeline, snp_id, abundancy, pool)
            # results['frequency'] = 
            all_results.append(results)
        except Exception as e:
            pipeline.print2('An exception happened:{}'.format(e))
        count += 1
        if count % 2 == 0:
            with open(os.path.join(pipeline.FILES['out_directory'], 'mixing_analysis_partial.json'), 'w') as f:
                json.dump(all_results, f, indent=3)

    try:
        os.remove(os.path.join(pipeline.FILES['out_directory'], 'mixing_analysis_partial.json'))
    except Exception as e:
        pass

    pipeline.print2('Mixing analysis complete. Now saving.')
    with open(os.path.join(pipeline.FILES['out_directory'], 'mixing_analysis.json'), 'w') as f:
        json.dump(all_results, f, indent=3)
    
    pool.close()
    pool.join()
    
    pipeline.print2('Mixing analysis saving complete.')


def _snp_mixing_analysis(pipeline, snp_id, abundancy, pool=None):
    """
    Does the analysis for a single snp
    """
    if pool is None:
        private_pool = True
        pool = multiprocessing.Pool(CPU_COUNT,maxtasksperchild=1)
    else:
        private_pool = False

    tumor = pipeline.tumor
    sampler = pipeline.sampler

    # 1) find cells that contain snp_id: "special" cells
    pipeline.print2('Searching for cells with the snp.')
    genotype_idx, cell_ids, cells = tumor.cells_with_snp(snp_id, pool)
    # genotype_idx: the genotypes that contain the snp_id in question
    # cell_ids: a list of all cells that have one of the genotypes in genotype_idx
    # cells: the cell object

    # if there are too many, montecarlo this
    pipeline.print2('Found ' +str(len(cell_ids)) +' cells for SNP: '+str(snp_id))

    """
    In [12]: list_to_sample = [MyObject()] * 5000000
    In [13]: %timeit random.sample(list_to_sample, 100)
    10000 loops, best of 3: 30.5 \mus per loop
    In [15]: %yimeit [list_to_sample[i] for i in np.random.randint(0, 5000000, 100)]
    100000 loops, best of 3: 17.3 \mus per loop

    Using numpy to generate indicies for sampling is much faster than using random.sample
    """
    if len(cell_ids) > 100:
        pipeline.print2('Monte Carlo Sampling')
        cell_ids = [cell_ids[i] for i in np.random.choice(len(cell_ids), 100)]
        # cell_ids = random.sample(cell_ids, 100)
    pipeline.print2('Sampled {} cells. Now performing mixing analysis'.format(len(cell_ids)))

    # obtain the "original_id" of the genotypes we have
    special_original_genotype_ids = set([g.original_id for g in tumor.get_genotypes(genotype_idx)])

    all_cells_COM = Statistics.centre_of_mass(cells).tolist()

    results = {
            'snp_id': snp_id,
            'COM': all_cells_COM,
            'frequency': abundancy
        }

    # special_proportion = []
    # COMS = []
    pipeline.print2('Preparing multiprocessing loop.')
    arguments_prepared = [(pipeline, cell_id, special_original_genotype_ids) for cell_id in cell_ids]
    pipeline.print2('Entering multiprocessing loop.')
    # for cell_id in cell_ids:
    #     _special_proportion = []
    #     # _COMS = []
    #     cell_position = tumor.cells[cell_id, 0:3]
    #     _, cell_positions_all, genotypes_all = sampler.sample_fixed_points(max(CLUSTER_SIZES), centre=cell_position)
    #     # 2) go through each "special cell" and look at small radius around it
    #     for n in CLUSTER_SIZES:
    #         # pick the first n cell positions:
    #         cell_positions = cell_positions_all[:n, :]
    #         sample_original_genotype_ids = set(map(lambda g : g.original_id, genotypes_all[:n]))

    #         special_genotypes_in_sample = sample_original_genotype_ids.intersection(special_original_genotype_ids)
    #         # print('specials in sample:',special_genotypes_in_sample)
    #         # print('all genotypes in samples',sample_original_genotype_ids)
    #         # print('specials for this snp', special_original_genotype_ids)
    #         assert len(special_genotypes_in_sample) >= 1, 'the sample must have at LEAST the original special cell'
    #         # 3) what proportion of it contains special cells? (can be done via genotype cross check)
    #         special_proportion_measured = len(special_genotypes_in_sample) / float(len(sample_original_genotype_ids))
            
    #         _special_proportion.append(special_proportion_measured)
        
    special_proportion = pool.map(_individual_cell_snp_mixing_analysis, arguments_prepared)
    pipeline.print2('Multiprocessing loop done.')
    
    special_proportion = np.array(special_proportion)
    std_special_proportion = special_proportion.std(axis=0).tolist()
    special_proportion = special_proportion.mean(axis=0).tolist()

    results['sample_sizes'] = CLUSTER_SIZES
    # results['COMs'] = COMS
    results['special_proportion'] = special_proportion
    results['std_special_proportion'] = std_special_proportion

    if private_pool:
        # only close these processes if we created them
        pool.close()
        pool.join()

    return results

def _individual_cell_snp_mixing_analysis(arguments):
    """
    Use for multiprocessing
    """
    pipeline, cell_id, special_original_genotype_ids = arguments
    # pipeline.print2('')
    sampler = pipeline.sampler
    tumor = pipeline.tumor

    _special_proportion = []
    # _COMS = []
    cell_position = tumor.cells[cell_id, 0:3]
    _, cell_positions_all, genotypes_all = sampler.sample_fixed_points(max(CLUSTER_SIZES), centre=cell_position)

    # pipeline.print2('CellID {}. Obtained a neighbourhood of size {}'.format(cell_id, len(cell_positions_all)))
    container_indicators = np.zeros(max(CLUSTER_SIZES), dtype=np.bool8)

    # this is O(max(CLUSTER_SIZES)) since set lookup is O(1)
    for i, genotype in enumerate(genotypes_all):
        container_indicators[i] = int(genotype.original_id in special_original_genotype_ids)
    # pipeline.print2('CellID {}. filled in indicators'.format(cell_id))
    # this is O(max(CLUSTER_SIZES))
    special_counts = np.cumsum(container_indicators)
    # pipeline.print2('CellID {}. cum sum calculated'.format(cell_id))
    assert special_counts[0] == 1, 'the first cell must be special by construction'
    for n in CLUSTER_SIZES:
        special_proportion_measured = special_counts[n-1] / float(n)
        _special_proportion.append(special_proportion_measured)

    # # 2) go through each "special cell" and look at small radius around it
    # for n in CLUSTER_SIZES:
    #     # pick the first n cell positions:
    #     cell_positions = cell_positions_all[:n, :]
    #     sample_original_genotype_ids = set(map(lambda g : g.original_id, genotypes_all[:n]))

    #     special_genotypes_in_sample = sample_original_genotype_ids.intersection(special_original_genotype_ids)
    #     # print('specials in sample:',special_genotypes_in_sample)
    #     # print('all genotypes in samples',sample_original_genotype_ids)
    #     # print('specials for this snp', special_original_genotype_ids)
    #     assert len(special_genotypes_in_sample) >= 1, 'the sample must have at LEAST the original special cell'
    #     # 3) what proportion of it contains special cells? (can be done via genotype cross check)
    #     special_proportion_measured = len(special_genotypes_in_sample) / float(len(sample_original_genotype_ids))
        
    
    return _special_proportion

