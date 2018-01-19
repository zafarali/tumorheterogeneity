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

EPS = np.finfo(float).eps
MAX_CELLS_RANGE = 50
INTERVAL = 5
LARGE_CLUSTERS = [100, 1000, 10000, 20000]
SMALL_CLUSTERS = range(2, MAX_CELLS_RANGE, INTERVAL)
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

    # 2) sample snps according to frequency
    mutation_idx = np.random.choice(mutations['Num'], size=50, replace=False, p=mutations['sampling_prob'])
    mutations_to_analyze = mutations.iloc[mutation_idx]

    # 3) perform analysis for each of those snps

    all_results = []
    for row in mutations_to_analyze.iterrows():
        data = row[1]
        abundancy, snp_id = data.abundancy, int(data.SNP)
        results = _snp_mixing_analysis(pipeline, snp_id)
        results['frequency'] = abundancy
        all_results.append(results)

    pipeline.print2('Mixing analysis complete. Now saving.')
    with open(os.path.join(pipeline.FILES['out_directory'], 'mixing_analysis.json'), 'w') as f:
        json.dump(all_results, f, indent=3)

    pipeline.print2('Mixing analysis saving complete.')


def _snp_mixing_analysis(pipeline, snp_id):
    """
    Does the analysis for a single snp
    """

    tumor = pipeline.tumor
    sampler = pipeline.sampler

    # 1) find cells that contain snp_id: "special" cells
    genotype_idx, cell_ids, cells = tumor.cells_with_snp(snp_id)
    # genotype_idx: the genotypes that contain the snp_id in question
    # cell_ids: a list of all cells that have one of the genotypes in genotype_idx
    # cells: the cell object

    # if there are too many, montecarlo this
    if len(cell_ids) > 100:
        pipeline.print2('Monte Carlo Sampling')
        cell_ids = random.sample(cell_ids, 100)
    
    # obtain the "original_id" of the genotypes we have
    special_original_genotype_ids = set([g.original_id for g in tumor.get_genotypes(genotype_idx)])

    all_cells_COM = Statistics.centre_of_mass(cells).tolist()

    results = {
            'snp_id': snp_id,
            'COM': all_cells_COM
        }

    special_proportion = []
    COMS = []

    pipeline.print2('Performing mixing analysis on ' +str(len(cell_ids)) +' cells for SNP: '+str(snp_id))
    for cell_id in cell_ids:
        _special_proportion = []
        _COMS = []
        cell_position = tumor.cells[cell_id, 0:3]
        _, cell_positions_all, genotypes_all = sampler.sample_fixed_points(max(CLUSTER_SIZES), centre=cell_position)
        # 2) go through each "special cell" and look at small radius around it
        for n in CLUSTER_SIZES:
            # pick the first n cell positions:
            cell_positions = cell_positions_all[:n, :]
            sample_original_genotype_ids = set(map(lambda g : g.original_id, genotypes_all[:n]))

            special_genotypes_in_sample = sample_original_genotype_ids.intersection(special_original_genotype_ids)
            # print('specials in sample:',special_genotypes_in_sample)
            # print('all genotypes in samples',sample_original_genotype_ids)
            # print('specials for this snp', special_original_genotype_ids)
            assert len(special_genotypes_in_sample) >= 1, 'the sample must have at LEAST the original special cell'
            # 3) what proportion of it contains special cells? (can be done via genotype cross check)
            special_proportion_measured = len(special_genotypes_in_sample) / float(len(sample_original_genotype_ids))
            
            COM = Statistics.centre_of_mass(cell_positions).tolist()
            _COMS.append(COM)
            _special_proportion.append(special_proportion_measured)

        special_proportion.append(_special_proportion)
        COMS.append(_COMS)

    special_proportion = np.array(special_proportion).mean(axis=0).tolist()
    COMS = np.array(COMS).mean(axis=0).tolist()

    results['sample_sizes'] = CLUSTER_SIZES
    results['COMs'] = COMS
    results['special_proportion'] = special_proportion

    return results





MIXING_MODULES = [prepare_tumor_for_mixing_analysis, perform_mixing_analysis]